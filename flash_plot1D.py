import numpy as np
import matplotlib.pyplot as plt
import yt
import matplotlib.colors as colors
import os
import glob
import quantity_dict
from make_hydro_movie import *
from format_hydro_files import *

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.size'] = 14

folderpath = "/Users/aodhan/Documents/HYDRO/magnetised_shocks/unmag_shock/"

filename = "unmag_shock_1_hdf5_plt_cnt_"

# [sim_vacuumHeight
centre = [4900] # scale the axis to put the centre of the spot to (0)
resolution = [500] # resolution of the fixed resolution buffer
n_crit = 1.01E21 # for normalising the densirt
folders = 'unmag_shock_1/'
files = range(0, 4)
plot_B = False
plot_bdry = False
plot_mass = True
eV_temp = True

print(dict)
if eV_temp:
    temp_factor=8.62E-5
else:
    temp_factor=1

def get_quantity(frb, quantity):
    data = (frb[quantity].d)
    data = data[~(data == 0).all(1)]
    return data

def central_lineout(x, y, quantites, labels, folderpath, filenumber, time, dictionary, save):
    fig, ax = plt.subplots()

    x_values = np.linspace(x[0], x[1], len(quantites[0]))
    counter = 0
    for quantity in quantites:
        ax.plot(x_values, quantity, label=labels[counter])
        counter+=1
    ax.set_title(time + 'ps')
    ax.set_yscale('log')
    ax.set_xlabel('x [$\mu$m]')
    ax.set_ylabel(dictionary["ylabel"])
    ax.set_ylim(dictionary["min"], dictionary["max"])
    ax.set_title(time + 'ps')
    fig.tight_layout()
    plt.legend(frameon=False)
    plt.savefig(folderpath + str(filenumber) + save, dpi=200)
    plt.close()

counter=0
for filenumber in files:

    hydro_file = format_file_number(filenumber)
    file = folderpath + folders + filename + hydro_file
    ds = yt.load(file)
    #cartesian
    slc = ds.slice('z', 0.5)
    print("Length unit: ", ds.length_unit)
    #cylindrical
    #slc = ds.slice(2, 0)

    time = str(round(ds.current_time.d*1E12, 2))
    if counter == 0:
        leftedge = ds.domain_left_edge.d*1E4
        rightedge = ds.domain_right_edge.d*1E4

        y = leftedge - centre, rightedge - centre
        
        target_centre = [rightedge[0] / 2E4, rightedge[1] / 2E4, 0] #centre of the simulation in code units (cm)
    # [x width in cm, y in cm]
    frb = slc.to_frb(width=(0.3, 'code_length'), resolution=resolution, height=(0.5, 'code_length'), center=target_centre)
    bounds = np.array(frb.bounds)
    x = (bounds[2] * 1E4) - centre[0], (bounds[3]*1E4) - centre[0]
    y = (bounds[0] * 1E4) - centre[1], (bounds[1]*1E4) - centre[1]
    
    density=get_quantity(frb, quantity='dens')
    YE=get_quantity(frb, quantity='ye')
    SUMY=get_quantity(frb, quantity='sumy')
    electron_temp=get_quantity(frb, quantity='tele')*temp_factor
    ion_temp=get_quantity(frb, quantity='tion')*temp_factor
    rad_temp=get_quantity(frb, quantity='trad')*temp_factor

    if plot_B:
        Bz = get_quantity(frb, quantity='magz')
    if plot_bdry:
        bdry = get_quantity(frb, quantity='bdry')
    
    electron_density = 6.023E23*YE*density
    # ionisation = YE/SUMY
    ion_density = 6.023E23*SUMY*density

    # 1d lineouts
    central_lineout(x, y, [electron_density, ion_density], ['e$^-$', 'ion'], folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_1D["density"], save='1Ddens')
    central_lineout(x, y, [electron_temp, ion_temp, rad_temp], ['e$^-$', 'ion', 'rad'], folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_1D["temp_eV"], save='1Dtemp')

    if plot_B:
        central_lineout(x, y, Bz, ['B$_z$'], folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_1D["magz"], save='Bz')

    if counter == len(files)-1:
        make_movie(folderpath + folders, files, '1Dtemp')
        make_movie(folderpath + folders, files, '1Ddens')
        if plot_B:
            make_movie(folderpath + folders, files, 'magz')
        files = glob.glob(folderpath + folders + '/*.png')
        for f in files:
            os.remove(f)

    counter +=1




