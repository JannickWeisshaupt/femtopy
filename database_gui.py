from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
import tkinter as tk
import matplotlib as mpl
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

class FixedDict(object):
    def __init__(self, dictionary):
        self._dictionary = dictionary

    def __setitem__(self, key, item):
        if key not in self._dictionary:
            raise KeyError("The key {} is not defined.".format(key))
        self._dictionary[key] = item

    def __getitem__(self, key):
        return self._dictionary[key]



class EmbeddedFigure(tk.Frame):
    def __init__(self, master, *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.master = master

        self.options_dict = FixedDict({'logarithmic': [False, False],'xlim': None,'ylim': None,'eV':False,'k':True})

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

    def update_figure(self, refractive_index_list=None):

        if refractive_index_list is None:
            if self.last_data is None:
                return
            else:
                refractive_index_list = self.last_data
        else:
            self.last_data = refractive_index_list

        self.subplot1.cla()
        self.subplot2.cla()

        x_min_list = []
        x_max_list = []
        y_min_list = []
        y_max_list = []

        color_cycle = mpl.rcParams['axes.prop_cycle']

        for refractive_index_handler,plot_color in zip(refractive_index_list,color_cycle):
            wav = refractive_index_handler.wavelength
            n = refractive_index_handler.n
            k = refractive_index_handler.k
            gvd  = refractive_index_handler.gvd

            if self.options_dict['eV']:
                x_plot = 1.240/wav
                self.subplot2.set_xlabel(r'Photon Energy [eV]')
            else:
                x_plot = wav
                self.subplot2.set_xlabel(r'Wavelength [$\mu$m]')

            x_min_list.append(x_plot.min())
            x_max_list.append(x_plot.max())
            if k is None or not self.options_dict['k']:
                y_min_list.append(n.min())
                y_max_list.append(n.max())
            else:
                y_min_list.append(min([n.min(),k.min()]))
                y_max_list.append(min([n.max(),k.max()]))


            if self.options_dict['xlim'] is not None:
                xmin = self.options_dict['xlim'][0]
                xmax = self.options_dict['xlim'][1]
                mask = (x_plot>=xmin) & (x_plot<=xmax)
                mask = np.roll(mask, 1) | mask | np.roll(mask, -1)
                n = n[mask]
                gvd = gvd[mask]
                if k is not None:
                    k = k[mask]
                wav = wav[mask]
                x_plot = x_plot[mask]

            plot_color=plot_color['color']
            self.subplot1.plot(x_plot,n,color=plot_color,label=refractive_index_handler.name)
            if k is not None and self.options_dict['k']:
                self.subplot1.plot(x_plot,k,'--',color=plot_color)

            self.subplot2.plot(x_plot,gvd,color=plot_color)

        if self.options_dict['xlim'] is None:
            self.subplot1.set_xlim(min(x_min_list), max(x_max_list))
        else:
            self.subplot1.set_xlim(self.options_dict['xlim'][0], self.options_dict['xlim'][1])


        if self.options_dict['ylim'] is not None:
            self.subplot1.set_ylim(self.options_dict['ylim'][0], self.options_dict['ylim'][1])

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

        if len(refractive_index_list)>1:
            self.subplot1.legend()


        self.subplot2.set_ylabel(r'GVD [fs$^2$/mm]')
        plt.pause(0.001)
        self.canvas.draw()


class RefrectiveIndexHandler:
    def __init__(self):
        self.n = None
        self.k = None
        self.wavelength = None
        self.gvd = None
        self.name = None
        self.c0 = 3e8

    def get_data_from_id(self,id):
        info = db.search_id(id)
        self.name = info['book']+' '+info['page']
        var1 = db.get_material_n_numpy(id)
        self.wavelength = var1[:,0]
        self.n = var1[:,1]
        var2 = db.get_material_k_numpy(id)
        if var2 is not None:
            self.k = var2[:,1]
        else:
            self.k = None
        wav2, gvd_res = self.calculate_gvd(self.wavelength, self.n)
        self.gvd = gvd_res

    def get_data_from_info(self,data):
        id = int(data[0])
        self.get_data_from_id(id)

    def derivative(self,n, order):
        l = n[:, 0]
        ref = n[:, 1]

        for i in range(order):
            ref = (ref - np.roll(ref, 1)) / (l - np.roll(l, 1))
        return np.array([l[order:], ref[order:]]).T

    def second_der(self,l, ref):
        var1 = np.roll(ref, -1) - 2 * ref + np.roll(ref, 1)
        var2 = l - np.roll(l, 1)
        res = var1 / var2 ** 2
        res[0] = np.nan
        res[-1] = np.nan
        return np.array([l, res]).T

    def calculate_gvd(self, l,ref):
        der_sol = self.second_der(l, ref)
        lambda0 = der_sol[:, 0]
        der_1 = der_sol[:, 1]
        res = lambda0 ** 3 / (2 * np.pi * self.c0 ** 2) * der_1
        res = res * 1e21
        return lambda0, res


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

        def delayed_click(event):
            master.after(10,lambda: self.click_event(event))

        self.result_tree.bind("<Button-1>", delayed_click)

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

    def click_event(self, event):
        selection_list = self.result_tree.selection()
        self.master.list_of_refractive_index = []
        for i,item_1 in enumerate(selection_list):
            res = self.result_tree.item(item_1,"values")
            if res == '':
                continue
            ref1 = RefrectiveIndexHandler()
            ref1.get_data_from_id(res[6])
            self.master.list_of_refractive_index.append(ref1)
            if i>20:
                tk.messagebox.showerror('Too many selected','You have more than twenty selected which is imho too much')
                break

        if len(self.master.list_of_refractive_index) > 0:
            q1.put(['update figure'])



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

        self.frame0 = tk.Frame(self)
        self.frame0.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.search_field = LabelWithEntry(self.frame0,'Search')
        self.search_field.pack(side=tk.LEFT)

        self.search_field.bind('<Key>',handle_update_event)

        self.exact_var = tk.IntVar()
        self.exact_var.set(False)
        self.exact_checkbutton = ttk.Checkbutton(self.frame0, text="Exact", variable=self.exact_var, takefocus=False,
                                                 command=handle_update_event)

        self.exact_checkbutton.pack(side=tk.LEFT, padx=5)

        self.search_options = {'Material': 'book','Data from':'page'}
        self.search_combo = NewCBox(self.frame0,self.search_options,current='Material')
        self.search_combo.bind("<<ComboboxSelected>>", handle_update_event)
        self.search_combo.pack(side=tk.LEFT, padx=5)

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
        def delayed_click(event):
            master.after(10,lambda: self.click_event(event))

        self.result_tree.bind("<Button-1>", delayed_click)



        self.result_tree.configure(yscrollcommand=self.scrollbar.set)

        self.search_res = None
        self.selected_res = None

    def search(self):
        exact_bool = self.exact_var.get()
        search_term = self.search_field.get()
        if search_term == '':
            self.result_tree.delete(*self.result_tree.get_children())
            return


        search_area = self.search_combo.value()
        if not exact_bool:
            sql_query = "SELECT * FROM pages WHERE "+search_area+" like '%"+search_term+"%'"
        else:
            sql_query = "SELECT * FROM pages WHERE "+search_area+" like '" + search_term + "'"
        self.search_res = db.search_custom(sql_query)
        # self.search_res = db.search_pages(search_term,exact=True)
        self.update_tree()

    def update_tree(self):
        self.result_tree.delete(*self.result_tree.get_children())
        for i,res in enumerate(self.search_res):
            self.result_tree.insert('', 'end', text=str(i),values=(res[2],res[3],res[7],res[8],('No','Yes')[res[6]],res[9]))

    def click_event(self, event):
        selection_list = self.result_tree.selection()
        self.master.list_of_refractive_index = []
        for item_1 in selection_list:
            item_index = int(self.result_tree.item(item_1, "text"))
            res = self.search_res[item_index]
            self.selected_res = res
            ref1 = RefrectiveIndexHandler()
            ref1.get_data_from_id(res[0])
            self.master.list_of_refractive_index.append(ref1)
        if len(selection_list) > 0:
            q1.put(['update figure'])


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

        self.k_var = tk.IntVar()
        self.k_var.set(self.option_dict['k'])
        self.k_checkbutton = ttk.Checkbutton(self, text="Show k", variable=self.k_var , takefocus=False,
                                                 command=self.set_values_for_figure)
        self.k_checkbutton.pack(side=tk.LEFT, padx=self.padx)

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

        if self.k_var.get():
            self.option_dict['k'] = True
        else:
            self.option_dict['k'] = False

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

        self.list_of_refractive_index = None

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
                self.embedded_figure.update_figure(self.list_of_refractive_index)

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
