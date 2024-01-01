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

# folderpath = "/Users/AMcilvenny/Documents/Hydro/magnetised_shocks/unmag_domenico/"
folderpath = '/Users/AMcilvenny/Desktop/'

filename = "TAW_2015_CSW_CH_hdf5_plt_cnt_"

# [sim_vacuumHeight, probably (xmax-xmin)/2] in microns
centre = [4900, 0] # scale the axis to put the centre of the spot to (0, 0)
# centre = [0, 0]
resolution = [250, 250] # resolution of the fixed resolution buffer
n_crit = 1.01E21 # for normalising the densirt
# folders = 'unmag21/'
folders = 'Conor/'
files = range(50, 51)
plot_B = False
plot_bdry = False
plot_mass = False
eV_temp = True
normalise_density = False
time_normalise = 0 # 1002.5 but ignoring the final few ps because these aren't modelled in flash
xlim = -5000
ylim = 1500

if eV_temp:
    temp_factor=8.62E-5
else:
    temp_factor=1

def get_quantity(frb, quantity):
    data = (frb[quantity].d).T
    data = data[~(data == 0).all(1)]
    return data

def plot_2D(x, y, quantity, folderpath, filenumber, time, dictionary):
    fig, ax = plt.subplots()   
    cax = ax.imshow(quantity, cmap=dictionary["cmap"],
             extent=[x[0], x[1], y[0], y[1]], norm=colors.LogNorm(vmin=dictionary["vmin"], vmax=dictionary["vmax"]), interpolation='nearest', origin='lower', aspect='auto')
    # cax = ax.imshow(quantity, cmap=dictionary["cmap"],
    #          extent=[x[0], x[1], y[0], y[1]], vmin=dictionary["vmin"], vmax=dictionary["vmax"], interpolation='nearest', origin='lower', aspect='auto')
    # ax.set_xlabel('x [$\mu$m]')
    ax.set_title(time + 'ps')
    ax.set_ylabel('y [$\mu$m]')
    ax.set_xlim(xlim, 0)
    ax.set_ylim(-ylim, ylim)
    # ax.set_xlim(-100, 10)
    # ax.set_ylim(-30, 30)
    cbar = fig.colorbar(cax)
    cbar.ax.set_ylabel(dictionary["cbar_label"])
    fig.tight_layout()
    plt.savefig(folderpath + filenumber + dictionary["name"] + '.png', dpi=200)
    plt.close()

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
    ax.set_xlim(xlim, 0)
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

    time = str(round(ds.current_time.d*1E12 - time_normalise, 2)) 
    if counter == 0:
        leftedge = ds.domain_left_edge.d*1E4
        rightedge = ds.domain_right_edge.d*1E4

        x = leftedge[0] - centre[0], rightedge[0] - centre[0]
        y = leftedge[1] - centre[1], rightedge[1] - centre[1]
        
        target_centre = [(leftedge[0] + rightedge[0]) / 2E4, (leftedge[1] + rightedge[1]) / 2E4, 0] #centre of the simulation in code units (cm)
    # [x width in cm, y in cm]
    frb = slc.to_frb(width=(0.3, 'code_length'), resolution=resolution, height=(0.5, 'code_length'), center=target_centre)
    bounds = np.array(frb.bounds)
    x = (bounds[2] * 1E4) - centre[0], (bounds[3]*1E4) - centre[0]
    y = (bounds[0] * 1E4) - centre[1], (bounds[1]*1E4) - centre[1]
    
    # standard quantities
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
    
    if normalise_density == True:
        density_factor = n_crit # convert to critical density
    else:
        density_factor =  1  # keeps in cm3
    electron_density = 6.023E23*YE*density/density_factor
    # ionisation = YE/SUMY
    ion_density = 6.023E23*SUMY*density/density_factor

    # 2d colorplots
    plot_2D(x, y, electron_temp, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_2D["tele"])
    plot_2D(x, y, ion_temp, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_2D["tion"])
    plot_2D(x, y, rad_temp, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_2D["trad"])
    plot_2D(x, y, electron_density, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_2D["edens_cm3"])
    plot_2D(x, y, ion_density, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_2D["iondens_cm3"])
    # plot_2D(x, y, electron_density, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_2D["edens"])
    # plot_2D(x, y, ion_density, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_2D["iondens"])

    # 1d lineouts
    central_lineout(x, y, [electron_density, ion_density], ['e$^-$', 'ion'], folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_1D["density_cm3"], save='1Ddens')
    central_lineout(x, y, [electron_temp, ion_temp, rad_temp], ['e$^-$', 'ion', 'rad'], folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_1D["temp_eV"], save='1Dtemp')

    # if plot_B:
    #     plot_2D(x, y, Bz, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict_2D["magz"])

    if counter == len(files)-1:
        make_movie(folderpath + folders, files, 'iondens')
        make_movie(folderpath + folders, files, 'edens')
        make_movie(folderpath + folders, files, 'tele')
        make_movie(folderpath + folders, files, 'tion')
        make_movie(folderpath + folders, files, 'trad')
        make_movie(folderpath + folders, files, '1Dtemp')
        make_movie(folderpath + folders, files, '1Ddens')
        if plot_B:
            make_movie(folderpath + folders, files, 'magz')
        files = glob.glob(folderpath + folders + '/*.png')
        for f in files:
            os.remove(f)

    counter +=1




