from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
import tkinter as tk
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib import gridspec
import numpy as np
import threading
import queue
import os
import sys
import ctypes

plt.rc('font', **{'size': 14})

q1 = queue.Queue()

from refractivesqlite import dboperations as DB
import matplotlib.pyplot as plt
import numpy as np

dbpath = r"RefractiveIndex.db"

db = DB.Database(dbpath)



class NewCBox(ttk.Combobox):
    def __init__(self, master, dictionary, current=0, *args, **kw):
        ttk.Combobox.__init__(self, master, values=list(dictionary.keys()), state='readonly', *args, **kw)
        self.dictionary = dictionary
        self.set(current)

    def value(self):
        return self.dictionary[self.get()]


class LabelWithEntry(tk.Frame):
    def __init__(self, master, text, width_label=10, width_entry=20, *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.label = ttk.Label(self, text=text, width=width_label)
        self.label.grid(row=0, column=0)

        self.entry = ttk.Entry(self, width=width_entry)
        self.entry.grid(row=0, column=1)

    def get(self):
        return self.entry.get()

    def set(self, value, fmt="{0:6.4f}"):
        self.entry.delete(0, 'end')
        self.entry.insert(0, fmt.format(value))

    def bind(self, sequence=None, func=None, add=None):
        self.entry.bind(sequence=sequence, func=func, add=add)


class EmbeddedFigure(tk.Frame):
    def __init__(self, master, *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.master = master

        self.options_dict = {'logarithmic': [False, False],'xlim': None}

        self.f = plt.Figure(figsize=(10, 6))
        gs = gridspec.GridSpec(1, 1, height_ratios=[1])
        self.subplot1 = self.f.add_subplot(gs[0])
        self.subplot1.format_coord = lambda x, y: "x={0:1.2f}, y={1:1.1f}".format(x, y)



        colortuple = master.winfo_rgb(self.master.cget('bg'))
        color_rgb = [x / 16 ** 4 for x in colortuple]
        self.f.patch.set_facecolor(color_rgb)
        self.subplot1.set_ylim(0, 1)

        self.subplot1.set_xlabel(r'Wavelength [$\mu$m]')
        self.subplot1.set_ylabel(r'Refractive Index')

        self.f.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.f, self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.last_data = None

    def update_figure(self, data=None):

        if data is None:
            if self.last_data is None:
                return
            else:
                data = self.last_data
        else:
            self.last_data = data

        id = int(data[0])
        n = db.get_material_n_numpy(id)
        k = db.get_material_k_numpy(id)

        self.subplot1.cla()

        wav=n[:,0]

        self.subplot1.plot(n[:,0],n[:,1])
        if k is not None:
            self.subplot1.plot(k[:,0],k[:,1])

        self.subplot1.set_xlim(wav.min(),wav.max())
        self.subplot1.set_xlabel(r'Wavelength [$\mu$m]')
        self.subplot1.set_ylabel(r'Refractive Index')

        plt.pause(0.001)
        self.canvas.draw()


class SearchFrame(tk.Frame):
    def __init__(self, master, *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.master = master

        self.search_field = LabelWithEntry(self,'Search')
        self.search_field.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.search_field.bind('<Return>',lambda event: self.search())

        self.frame1 = tk.Frame(self)
        self.frame1.pack(side=tk.TOP, fill=tk.BOTH, expand=True)



        self.result_tree = ttk.Treeview(self.frame1, columns=('Index','Material','Data from','rangeMin', 'rangeMax','Extinction','Points'))

        self.scrollbar = ttk.Scrollbar(self.frame1, orient="vertical", command=self.result_tree.yview)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.result_tree.heading('#0', text='Index', anchor=tk.CENTER)
        self.result_tree.column("#0", minwidth=50, width=50, stretch=tk.NO)
        self.result_tree.heading('#1', text='Material', anchor=tk.CENTER)
        self.result_tree.column("#1", minwidth=50, width=200, stretch=tk.NO)
        self.result_tree.heading('#2', text='Data from', anchor=tk.CENTER)
        self.result_tree.column("#2", minwidth=50, width=200, stretch=tk.NO)
        self.result_tree.heading('#3', text='rangeMin', anchor=tk.CENTER)
        self.result_tree.column("#3", minwidth=50, width=130, stretch=tk.NO)
        self.result_tree.heading('#4', text='rangeMax', anchor=tk.CENTER)
        self.result_tree.column("#4", minwidth=50, width=130, stretch=tk.NO)
        self.result_tree.heading('#5', text='Extinction', anchor=tk.CENTER)
        self.result_tree.column("#5", minwidth=50, width=130, stretch=tk.NO)
        self.result_tree.heading('#6', text='Points', anchor=tk.CENTER)
        self.result_tree.column("#6", minwidth=50, width=130, stretch=tk.NO)

        self.result_tree.pack(side=tk.LEFT)
        self.result_tree.bind("<Double-1>", self.OnDoubleClick)



        self.result_tree.configure(yscrollcommand=self.scrollbar.set)

        self.search_res = None
        self.selected_res = None

    def search(self):
        exact_bool = False
        search_term = self.search_field.get()
        if not exact_bool:
            sql_query = "SELECT * FROM pages WHERE book like '%"+search_term+"%'"
        else:
            sql_query = "SELECT * FROM pages WHERE book like '" + search_term + "'"
        self.search_res = db.search_custom(sql_query)
        # self.search_res = db.search_pages(search_term,exact=True)
        self.update_tree()

    def update_tree(self):
        self.result_tree.delete(*self.result_tree.get_children())
        for i,res in enumerate(self.search_res):
            self.result_tree.insert('', 'end', text=str(i),values=(res[2],res[3],res[7],res[8],res[6],res[9]))

    def OnDoubleClick(self, event):
        item_1 = self.result_tree.selection()[0]
        item_index = int(self.result_tree.item(item_1, "text"))
        res = self.search_res[item_index]
        self.selected_res = res
        q1.put(['update figure',res])

class OptionWindow(tk.Toplevel):
    def __init__(self, master, *args, **kwargs):
        super().__init__(master, *args, **kwargs)

        self.master = master

        self.time_column_entry = LabelWithEntry(self, 'Column time:', width_label=20, width_entry=5)
        self.time_column_entry.set(0, fmt="{0:2d}")
        self.time_column_entry.pack()

        self.thz_column_entry = LabelWithEntry(self, 'Column Signal:', width_label=20, width_entry=5)
        self.thz_column_entry.set(0, fmt="{0:2d}")
        self.thz_column_entry.pack()

        self.ok_button = ttk.Button(self, text="OK", command=self.set_values, takefocus=False)
        self.ok_button.pack()

    def set_values(self):
        print(int(self.time_column_entry.get()))
        print(int(self.thz_column_entry.get()))



class FigureOptionFrame(tk.Frame):
    def __init__(self, master, option_dict, *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.master = master
        self.option_dict = option_dict

        self.padx = 2

        self.log_var_x = tk.IntVar()
        self.log_var_x.set(self.option_dict['logarithmic'][0])
        self.log_checkbutton_x = ttk.Checkbutton(self, text="x-Log", variable=self.log_var_x, takefocus=False,
                                                 command=self.set_values_for_figure)
        self.log_checkbutton_x.pack(side=tk.LEFT, padx=self.padx)

        self.log_var_y = tk.IntVar()
        self.log_var_y.set(self.option_dict['logarithmic'][1])
        self.log_checkbutton_y = ttk.Checkbutton(self, text="y-Log", variable=self.log_var_y, takefocus=False,
                                                 command=self.set_values_for_figure)
        self.log_checkbutton_y.pack(side=tk.LEFT, padx=self.padx)

        self.xmin_entry = LabelWithEntry(self,'fmin=', width_entry=5,width_label=5)
        self.xmin_entry.bind("<Return>",lambda event: self.set_values_for_figure())
        self.xmin_entry.pack(side=tk.LEFT, padx=self.padx)

        self.xmax_entry = LabelWithEntry(self,'fmax=', width_entry=5,width_label=5)
        self.xmax_entry.bind("<Return>",lambda event: self.set_values_for_figure())
        self.xmax_entry.pack(side=tk.LEFT, padx=self.padx)

    def set_values_for_figure(self):
        try:
            self.option_dict['xlim'] = [float(self.xmin_entry.get()), float(self.xmax_entry.get())]
        except ValueError:
            self.option_dict['xlim'] = None

        if self.log_var_x.get():
            self.option_dict['logarithmic'][0] = True
        else:
            self.option_dict['logarithmic'][0] = False

        if self.log_var_y.get():
            self.option_dict['logarithmic'][1] = True
        else:
            self.option_dict['logarithmic'][1] = False

        q1.put('update figure')


class Application(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.pack()
        self.master = master

        # master.resizable(width=False, height=False)
        master.geometry('{}x{}'.format(1000, 920))

        # master.wm_attributes("-topmost", 1)

        self.search_frame = SearchFrame(self)
        self.search_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.embedded_figure = EmbeddedFigure(self)
        self.embedded_figure.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.figure_option_frame = FigureOptionFrame(self, self.embedded_figure.options_dict)
        self.figure_option_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.option_window = None

        self.create_menubar()
        self.bind_global_commands()

        master.title('Example Gui')
        master.config(menu=self.menubar)
        master.protocol('WM_DELETE_WINDOW', self.quit_program)
        master.after(100, self.update_status)

    def update_status(self):
        root.after(100, self.update_status)
        if not q1.empty():
            q_el = q1.get()
            if q_el[0] == 'update figure':
                self.embedded_figure.update_figure(q_el[1])

    def create_menubar(self):
        self.menubar = tk.Menu(self, tearoff=0)
        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="File", menu=self.filemenu)

        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=self.quit_program, accelerator="Strq+q")

        self.optionmenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Options", menu=self.optionmenu)

        self.optionmenu.add_command(label='Example', command=self.open_option_window)

    def quit_program(self):
        self.master.destroy()

    def bind_global_commands(self):
        self.bind_all('<Control-q>', lambda event: self.quit_program())

    def open_option_window(self):
        if self.option_window is None or not self.option_window.winfo_exists():
            self.option_window = OptionWindow(self)


root = tk.Tk()
# root.iconbitmap(r'icon.ico')
app = Application(master=root)
root.mainloop()
