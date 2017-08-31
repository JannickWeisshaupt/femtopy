from cx_Freeze import setup, Executable, hooks
import os
import scipy
# use by python compile_setup.py build

os.environ['TCL_LIBRARY'] = r"C:\Users\weisshau\AppData\Local\Programs\Python\Python36-32\tcl\tcl8.6"
os.environ['TK_LIBRARY'] = r"C:\Users\weisshau\AppData\Local\Programs\Python\Python36-32\tcl\tk8.6"
# NOTE: you can include any other necessary external imports here aswell

# If you delete the whole folder you need to copy all .dll files from C:\Python34\Lib\site-packages\numpy\core
# into the folder of the executable (build/exe.win32-3.4)

includefiles = [r'C:\Users\weisshau\AppData\Local\Programs\Python\Python36-32\DLLs\tk86t.dll',r'C:\Users\weisshau\AppData\Local\Programs\Python\Python36-32\DLLs\tcl86t.dll',
                r"D:\own_python_modules\pyterminal\TerminalClass.py",r"D:\own_python_modules\pyterminal\terminal_frame.py","RefractiveIndex.db",'icon.ico']  # include any files here that you wish

scipy_path = os.path.dirname(scipy.__file__)
includefiles.append(scipy_path)

includes = ['numpy.core._methods','numpy.lib.format','numpy.matlib']
excludes = ['PyQt5','PyQt4','mpl_toolkits.basemap']
packages = ['pygments', 'pygments.lexers']

exe = Executable(
    # what to build
    script="database_gui.py",  # the name of your main python script goes here
    initScript=None,
    base='Win32GUI',  # if creating a GUI instead of a console app, type "Win32GUI"
    targetName="Refractive_Index.exe",  # this is the name of the executable file
    icon="./icon.ico"  # if you want to use an icon file, specify the file name here
)

setup(
    # the actual setup & the definition of other misc. info
    name="Refractive Index Database Viewer",  # program name
    version="0.0",
    description='',
    author="Jannick Weisshaupt",
    options={"build_exe": {"excludes": excludes, "packages": packages,
                           "include_files": includefiles, "includes": includes}},
    executables=[exe]
)