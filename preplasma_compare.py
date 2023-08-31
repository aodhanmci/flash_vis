import matplotlib.pyplot as plt
import numpy as np

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

folderpath = "/Users/AMcilvenny/Documents/Hydro/ip2_preplasma/"

filename = "laser_dim_hdf5_plt_cnt_"

# [sim_vacuumHeight, probably (xmax-xmin)/2] in microns
centre = [70, 100] # scale the axis to put the centre of the spot to (0, 0)
# centre = [0, 0]
resolution = [1000, 1000] # resolution of the fixed resolution buffer
n_crit = 1.75E21 # for normalising the densirt
folders = ['pre_4/', 'pre_6/']
laser_energy = ['5J', '17J']
plot_B = False
plot_bdry = False
plot_mass = True
eV_temp = True
normalise_density = True

if eV_temp:
    temp_factor=8.62E-5
else:
    temp_factor=1

def get_quantity(frb, quantity):
    data = (frb[quantity].d).T
    data = data[~(data == 0).all(1)]
    return data

def central_lineout(x, y, quantites, labels, folderpath, filenumber, time, dictionary, save):
    fig, ax = plt.subplots()
    if resolution[0] % 2 == 1:
        mid = int((resolution[0] + 1) / 2)
    else:
        mid = int(resolution[0] / 2)
    x_values = np.linspace(x[0], x[1], len(quantites[0]))
    counter = 0
    for quantity in quantites:
        data = quantity[mid, :]
        ax.plot(x_values, data, label=labels[counter])
        counter+=1
    ax.set_title(time + 'ps')
    # ax.set_xlim(-100, 5)
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
fig, ax = plt.subplots()
for folder in folders:

    hydro_file = format_file_number(20)
    file = folderpath + folder + filename + hydro_file

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

        x = leftedge[0] - centre[0], rightedge[0] - centre[0]
        y = leftedge[1] - centre[1], rightedge[1] - centre[1]
        
        target_centre = [(leftedge[0] + rightedge[0]) / 2E4, (leftedge[1] + rightedge[1]) / 2E4, 0] #centre of the simulation in code units (cm)
    # [x width in cm, y in cm]
    frb = slc.to_frb(width=(0.02, 'code_length'), resolution=resolution, height=(0.01, 'code_length'), center=target_centre)
    bounds = np.array(frb.bounds)
    x = (bounds[2] * 1E4) - centre[0], (bounds[3]*1E4) - centre[0]
    y = (bounds[0] * 1E4) - centre[1], (bounds[1]*1E4) - centre[1]
    
    # standard quantities
    density=get_quantity(frb, quantity='dens')
    YE=get_quantity(frb, quantity='ye')
    SUMY=get_quantity(frb, quantity='sumy')

    
    if normalise_density == True:
        density_factor = n_crit # convert to critical density
    else:
        density_factor =  1  # keeps in cm3
    electron_density = 6.023E23*YE*density/density_factor
    # ionisation = YE/SUMY
    # ion_density = 6.023E23*SUMY*density/density_factor


    
    if resolution[0] % 2 == 1:
        mid = int((resolution[0] + 1) / 2)
    else:
        mid = int(resolution[0] / 2)
    x_values = np.linspace(x[0], x[1], 1000)

    data = electron_density[mid, :]
    ax.plot(x_values, data, label=laser_energy[counter])
    counter+=1

ax.set_title(time + 'ps')
ax.set_yscale('log')
ax.set_xlabel('x [$\mu$m]')
ax.set_ylabel('n_e [n_c]')
ax.set_xlim(-70, 20)
ax.set_ylim(1e-2, 10)
ax.set_title(time + 'ps')
fig.tight_layout()
plt.legend(frameon=False)
plt.savefig('foam.png', dpi=200)
plt.close()





