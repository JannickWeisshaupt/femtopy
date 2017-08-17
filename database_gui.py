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
from refractivesqlite import dboperations as DB
import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', **{'size': 14})

q1 = queue.Queue()

dbpath = r"RefractiveIndex.db"

db = DB.Database(dbpath)

c0 = 3e8


def derivative(n,order):
    l = n[:,0]
    ref = n[:,1]

    for i in range(order):
        ref = (ref-np.roll(ref,1))/(l-np.roll(l,1))
    return np.array([l[order:],ref[order:]]).T

def second_der(n):
    l = n[:,0]
    ref = n[:,1]

    var1 = np.roll(ref,-1) - 2*ref+np.roll(ref,1)
    var2 = l-np.roll(l,1)
    res = var1/var2**2
    return np.array([l[1:-1],res[1:-1]]).T


def gvd(data):
    der_sol = second_der(data)
    lambda0 = der_sol[:,0]
    der_1 = der_sol[:,1]
    res = lambda0**3/(2*np.pi*c0**2)*der_1
    res = res * 1e21
    return lambda0,res


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

        self.options_dict = {'logarithmic': [False, False],'xlim': None,'ylim': None,'eV':False,'cum':False}

        self.f = plt.Figure(figsize=(10, 6))
        gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
        self.subplot1 = self.f.add_subplot(gs[0])
        self.subplot2 = self.f.add_subplot(gs[1],sharex=self.subplot1)
        self.subplot1.format_coord = self.formatter1
        self.subplot2.format_coord = lambda x, y: "lambda={0:1.2f}, GVD={1:1.0f}".format(x, y)


        colortuple = master.winfo_rgb(self.master.cget('bg'))
        color_rgb = [x / 16 ** 4 for x in colortuple]
        self.f.patch.set_facecolor(color_rgb)
        self.subplot1.set_ylim(0, 1)
        plt.setp(self.subplot1.get_xticklabels(), visible=False)

        self.subplot2.set_xlabel(r'Wavelength [$\mu$m]')
        self.subplot1.set_ylabel(r'Refractive Index')
        self.subplot2.set_ylim(0,1000)
        self.subplot2.set_ylabel(r'GVD [fs$^2$/mm]')

        self.f.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.f, self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.last_data = None

    def formatter1(self,x,y):
        return "lambda={0:1.2f}, n/k={1:1.3f}".format(x, y)

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
        self.subplot2.cla()

        wav=n[:,0]
        wav2,gvd_res = gvd(n)

        if self.options_dict['eV']:
            x_plot = 1.240/wav
            x_plot2 = 1.240/wav2
            self.subplot2.set_xlabel(r'Photon Energy [eV]')
        else:
            x_plot = wav
            x_plot2 = wav2
            self.subplot2.set_xlabel(r'Wavelength [$\mu$m]')



        if self.options_dict['xlim'] is None:
            self.subplot1.set_xlim(x_plot.min(), x_plot.max())
        else:
            self.subplot1.set_xlim(self.options_dict['xlim'][0], self.options_dict['xlim'][1])
            xmin = self.options_dict['xlim'][0]
            xmax = self.options_dict['xlim'][1]
            mask = (x_plot>=xmin) & (x_plot<=xmax)
            mask2 = (x_plot2>=xmin) & (x_plot2<=xmax)

            n = n[mask,:]
            gvd_res = gvd_res[mask2]
            if k is not None:
                k = k[mask,:]
            wav = wav[mask]
            x_plot = x_plot[mask]
            x_plot2 = x_plot2[mask2]



        if self.options_dict['ylim'] is not None:
            self.subplot1.set_ylim(self.options_dict['ylim'][0], self.options_dict['ylim'][1])



        self.subplot1.plot(x_plot,n[:,1])
        if k is not None:
            self.subplot1.plot(x_plot,k[:,1])

        plt.setp(self.subplot1.get_xticklabels(), visible=False)
        self.subplot1.set_ylabel(r'Refractive Index')

        if self.options_dict['logarithmic'][0]:
            self.subplot1.set_xscale('log')
        else:
            self.subplot1.set_xscale('linear')

        if self.options_dict['logarithmic'][1]:
            self.subplot1.set_yscale('log')
        else:
            self.subplot1.set_yscale('linear')

        self.subplot2.plot(x_plot2,gvd_res)
        self.subplot2.set_ylabel(r'GVD [fs$^2$/mm]')
        plt.pause(0.001)
        self.canvas.draw()


class WholeDatabaseFrame(tk.Toplevel):

    def __init__(self,master,*args,**kwargs):
        super().__init__(master, *args, **kwargs)
        self.master = master
        self.geometry('{}x{}'.format(1300, 600))
        self.frame1 = tk.Frame(self)
        self.frame1.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        def treeview_sort_column(tv, col, reverse):
            l = [(tv.set(k, col), k) for k in tv.get_children('')]
            l.sort(reverse=reverse)

            # rearrange items in sorted positions
            for index, (val, k) in enumerate(l):
                tv.move(k, '', index)

            # reverse sort next time
            tv.heading(col, command=lambda: \
                treeview_sort_column(tv, col, not reverse))

        self.result_tree = ttk.Treeview(self.frame1, columns=(
        'Index', 'Material', 'Data from', 'rangeMin', 'rangeMax', 'Extinction', 'Points','Db number'))

        self.scrollbar = ttk.Scrollbar(self.frame1, orient="vertical", command=self.result_tree.yview)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.BOTH)

        self.result_tree.column("#0", minwidth=50, width=150, stretch=tk.NO)
        self.result_tree.heading('#0', text='Index', anchor=tk.CENTER)
        self.result_tree.heading('#1', text='Material', anchor=tk.CENTER,
                                 command=lambda: treeview_sort_column(self.result_tree, '#1', False))
        self.result_tree.column("#1", minwidth=50, width=200, stretch=tk.NO)
        self.result_tree.heading('#2', text='Data from', anchor=tk.CENTER,
                                 command=lambda: treeview_sort_column(self.result_tree, '#2', False))
        self.result_tree.column("#2", minwidth=50, width=200, stretch=tk.NO)
        self.result_tree.heading('#3', text='rangeMin', anchor=tk.CENTER,
                                 command=lambda: treeview_sort_column(self.result_tree, '#3', False))
        self.result_tree.column("#3", minwidth=50, width=130, stretch=tk.NO)
        self.result_tree.heading('#4', text='rangeMax', anchor=tk.CENTER,
                                 command=lambda: treeview_sort_column(self.result_tree, '#4', False))
        self.result_tree.column("#4", minwidth=50, width=130, stretch=tk.NO)
        self.result_tree.heading('#5', text='Extinction', anchor=tk.CENTER,
                                 command=lambda: treeview_sort_column(self.result_tree, '#5', False))
        self.result_tree.column("#5", minwidth=50, width=130, stretch=tk.NO)
        self.result_tree.heading('#6', text='Points', anchor=tk.CENTER,
                                 command=lambda: treeview_sort_column(self.result_tree, '#6', False))
        self.result_tree.column("#6", minwidth=50, width=130, stretch=tk.NO)

        self.result_tree.heading('#7', text='Db number', anchor=tk.CENTER,
                                 command=lambda: treeview_sort_column(self.result_tree, '#6', False))
        self.result_tree.column("#7", minwidth=50, width=130, stretch=tk.NO)

        self.result_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.result_tree.configure(yscrollcommand=self.scrollbar.set)
        self.result_tree.bind("<Double-1>", self.OnDoubleClick)

        shelf_list = ['main','organic','glass','other','3d']
        for shelf in shelf_list:
            main = db.search_custom("SELECT * FROM pages WHERE shelf LIKE '"+shelf+"'")

            set_of_books = set()

            for el in main:
                set_of_books.add(el[2])

            x1 = self.result_tree.insert('', 'end', text=shelf,values=())
            for book in set_of_books:
                s_res = db.search_custom("SELECT * FROM pages WHERE shelf LIKE '"+shelf+"' and book like '"+book+"'")
                if len(s_res)>1:
                    y1 = self.result_tree.insert(x1, 'end', text=book,values=())
                    for i, res in enumerate(s_res):
                        self.result_tree.insert(y1, 'end', text='',
                                                    values=(res[2], res[3], res[7], res[8], ('No', 'Yes')[res[6]], res[9], res[0]))
                else:
                    res = s_res[0]
                    self.result_tree.insert(x1, 'end', text='', values=(res[2], res[3], res[7], res[8], ('No', 'Yes')[res[6]], res[8], res[0]))

    def OnDoubleClick(self, event):
        try:
            item_1 = self.result_tree.selection()[0]
            item_index = int(self.result_tree.item(item_1, "values")[6])
            res = db.search_id(item_index)
            q1.put(['update figure', list(res.values())])
        except Exception:
            pass

class SearchFrame(tk.Frame):
    def __init__(self, master, *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.master = master
        self._after_id = None

        def handle_update_event(*args):
            # cancel the old job
            if self._after_id is not None:
                master.after_cancel(self._after_id)
            # create a new job
            self._after_id = master.after(20, self.search)

        self.search_field = LabelWithEntry(self,'Search')
        self.search_field.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.search_field.bind('<Key>',handle_update_event)

        self.search_field.bind('<Return>',lambda event: self.search())



        self.frame1 = tk.Frame(self)
        self.frame1.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        def treeview_sort_column(tv, col, reverse):
            l = [(tv.set(k, col), k) for k in tv.get_children('')]
            l.sort(reverse=reverse)

            # rearrange items in sorted positions
            for index, (val, k) in enumerate(l):
                tv.move(k, '', index)

            # reverse sort next time
            tv.heading(col, command=lambda: \
                treeview_sort_column(tv, col, not reverse))

        self.result_tree = ttk.Treeview(self.frame1, columns=('Index','Material','Data from','rangeMin', 'rangeMax','Extinction','Points'))



        self.scrollbar = ttk.Scrollbar(self.frame1, orient="vertical", command=self.result_tree.yview)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.result_tree.column("#0", minwidth=50, width=50, stretch=tk.NO)
        self.result_tree.heading('#0', text='Index', anchor=tk.CENTER)
        self.result_tree.heading('#1', text='Material', anchor=tk.CENTER,command=lambda: treeview_sort_column(self.result_tree, '#1', False))
        self.result_tree.column("#1", minwidth=50, width=200, stretch=tk.NO)
        self.result_tree.heading('#2', text='Data from', anchor=tk.CENTER,command=lambda: treeview_sort_column(self.result_tree, '#2', False))
        self.result_tree.column("#2", minwidth=50, width=200, stretch=tk.NO)
        self.result_tree.heading('#3', text='rangeMin', anchor=tk.CENTER,command=lambda: treeview_sort_column(self.result_tree, '#3', False))
        self.result_tree.column("#3", minwidth=50, width=130, stretch=tk.NO)
        self.result_tree.heading('#4', text='rangeMax', anchor=tk.CENTER,command=lambda: treeview_sort_column(self.result_tree, '#4', False))
        self.result_tree.column("#4", minwidth=50, width=130, stretch=tk.NO)
        self.result_tree.heading('#5', text='Extinction', anchor=tk.CENTER,command=lambda: treeview_sort_column(self.result_tree, '#5', False))
        self.result_tree.column("#5", minwidth=50, width=130, stretch=tk.NO)
        self.result_tree.heading('#6', text='Points', anchor=tk.CENTER,command=lambda: treeview_sort_column(self.result_tree, '#6', False))
        self.result_tree.column("#6", minwidth=50, width=130, stretch=tk.NO)

        self.result_tree.pack(side=tk.LEFT)
        self.result_tree.bind("<Double-1>", self.OnDoubleClick)



        self.result_tree.configure(yscrollcommand=self.scrollbar.set)

        self.search_res = None
        self.selected_res = None

    def search(self):
        exact_bool = False
        search_term = self.search_field.get()
        if search_term == '':
            self.result_tree.delete(*self.result_tree.get_children())
            return

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
            self.result_tree.insert('', 'end', text=str(i),values=(res[2],res[3],res[7],res[8],('No','Yes')[res[6]],res[9]))

    def OnDoubleClick(self, event):
        item_1 = self.result_tree.selection()[0]
        item_index = int(self.result_tree.item(item_1, "text"))
        res = self.search_res[item_index]
        self.selected_res = res
        q1.put(['update figure',res])


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

        self.ev_var = tk.IntVar()
        self.ev_var.set(self.option_dict['eV'])
        self.ev_checkbutton = ttk.Checkbutton(self, text="eV", variable=self.ev_var, takefocus=False,
                                                 command=self.set_values_for_figure)
        self.ev_checkbutton.pack(side=tk.LEFT, padx=self.padx)

        self.xmin_entry = LabelWithEntry(self,'xmin=', width_entry=5,width_label=5)
        self.xmin_entry.bind("<Return>",lambda event: self.set_values_for_figure())
        self.xmin_entry.pack(side=tk.LEFT, padx=self.padx)

        self.xmax_entry = LabelWithEntry(self,'xmax=', width_entry=5,width_label=5)
        self.xmax_entry.bind("<Return>",lambda event: self.set_values_for_figure())
        self.xmax_entry.pack(side=tk.LEFT, padx=self.padx)

        self.ymin_entry = LabelWithEntry(self,'ymin=', width_entry=5,width_label=5)
        self.ymin_entry.bind("<Return>",lambda event: self.set_values_for_figure())
        self.ymin_entry.pack(side=tk.LEFT, padx=self.padx)

        self.ymax_entry = LabelWithEntry(self,'ymax=', width_entry=5,width_label=5)
        self.ymax_entry.bind("<Return>",lambda event: self.set_values_for_figure())
        self.ymax_entry.pack(side=tk.LEFT, padx=self.padx)

        self.cum_var = tk.IntVar()
        self.cum_var.set(self.option_dict['cum'])
        self.cum_checkbutton = ttk.Checkbutton(self, text="Cumulative", variable=self.cum_var , takefocus=False,
                                                 command=self.set_values_for_figure)
        self.cum_checkbutton.pack(side=tk.LEFT, padx=self.padx)

    def set_values_for_figure(self):
        try:
            self.option_dict['xlim'] = [float(self.xmin_entry.get()), float(self.xmax_entry.get())]
        except ValueError:
            self.option_dict['xlim'] = None

        try:
            self.option_dict['ylim'] = [float(self.ymin_entry.get()), float(self.ymax_entry.get())]
        except ValueError:
            self.option_dict['ylim'] = None

        if self.log_var_x.get():
            self.option_dict['logarithmic'][0] = True
        else:
            self.option_dict['logarithmic'][0] = False

        if self.log_var_y.get():
            self.option_dict['logarithmic'][1] = True
        else:
            self.option_dict['logarithmic'][1] = False

        if self.ev_var.get():
            self.option_dict['eV'] = True
        else:
            self.option_dict['eV'] = False

        if self.cum_var.get():
            self.option_dict['cum'] = True
        else:
            self.option_dict['cum'] = False

        q1.put(('update figure',None))


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

        self.whole_database_window = None

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

        self.datamenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Database", menu=self.datamenu)

        self.datamenu.add_command(label='Whole Database', command=self.open_whole_database_window)

    def quit_program(self):
        self.master.destroy()

    def bind_global_commands(self):
        self.bind_all('<Control-q>', lambda event: self.quit_program())

    def open_whole_database_window(self):
        if self.whole_database_window is None or not self.whole_database_window.winfo_exists():
            self.whole_database_window = WholeDatabaseFrame(self)


root = tk.Tk()
# root.iconbitmap(r'icon.ico')
app = Application(master=root)
root.mainloop()
