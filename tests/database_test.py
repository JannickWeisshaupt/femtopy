from refractivesqlite import dboperations as DB
import matplotlib.pyplot as plt
import numpy as np

dbpath = r"D:\own_python_modules\RefractiveIndex\RefractiveIndex.db"
ymlpath = r"D:\own_python_modules\RefractiveIndex"

db = DB.Database(dbpath)
db.create_database_from_folder(ymlpath,interpolation_points=500)

sres = db.search_pages("Si",exact=False)
n = db.get_material_n_numpy(61)
k = db.get_material_k_numpy(61)


def derivative(n,order):
    l = n[:,0]
    ref = n[:,1]

    for i in range(order):
        ref = (ref-np.roll(ref,1))/(l-np.roll(l,1))
    return np.array([l[order:],ref[order:]]).T

