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

# [sim_vacuumHeight, probably (xmax-xmin)/2] in microns
centre = [4900, 1500] # scale the axis to put the centre of the spot to (0, 0)
# centre = [0, 0]
resolution = [500, 500] # resolution of the fixed resolution buffer
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
    data = (frb[quantity].d).T
    data = data[~(data == 0).all(1)]
    return data

def plot_2D(x, y, quantity, folderpath, filenumber, time, dictionary):
    fig, ax = plt.subplots()
    cax = ax.imshow(quantity, cmap=dictionary["cmap"],
             extent=[x[0], x[1], y[0], y[1]], norm=colors.LogNorm(vmin=dictionary["vmin"], vmax=dictionary["vmax"]), interpolation='nearest', origin='lower', aspect='auto')
    ax.set_xlabel('x [$\mu$m]')
    ax.set_title(time + 'ps')
    ax.set_ylabel('y [$\mu$m]')
    cbar = fig.colorbar(cax)
    cbar.ax.set_ylabel(dictionary["cbar_label"])
    fig.tight_layout()
    plt.savefig(folderpath + filenumber + dictionary["name"] + '.png', dpi=200)
    plt.close()

def central_lineout_temp(x, y, e_temp, ion_temp, rad_temp, folderpath, filenumber, time, figure, axes):
    if resolution[0] % 2 == 1:
        mid = int((resolution[0] + 1) / 2)
    else:
        mid = int(resolution[0] / 2)
    # y_values = np.linspace(y[0], y[1], len(e_temp[:, 0]))
    # zero_line = np.argmin(np.abs(y_values))
    e_temp = e_temp[mid, :]
    ion_temp = ion_temp[mid, :]
    rad_temp = rad_temp[mid, :]
    x_values = np.linspace(x[0], x[1], len(e_temp))

    axes.plot(x_values, e_temp, label='e$^-$')
    axes.plot(x_values, ion_temp, label='ion')
    axes.plot(x_values, rad_temp, label='rad')
    axes.set_yscale('log')
    axes.set_xlabel('x [$\mu$m]')
    axes.set_ylabel('T [eV]')
    axes.legend()
    axes.set_title(time + 'ps')
    # axes.set_xlim(-50, 5)
    axes.set_ylim(1E-2, 1E4)
    # axes.set_xlim(-900, 100)
    axes.set_title(time + 'ns')
    plt.legend(frameon=False)
    plt.savefig(folderpath + str(filenumber) + '_centraltemp.png', dpi=200)
    plt.close()
    return x_values, e_temp, ion_temp


def central_lineout(x, y, electron_density, ion_density, folderpath, filenumber, time, figure, axes):
    if resolution[0] % 2 == 1:
        mid = int((resolution[0] + 1) / 2)
    else:
        mid = int(resolution[0] / 2)

    electron_density = electron_density[mid, :]
    ion_density = ion_density[mid, :]
    x_values = np.linspace(x[0], x[1], len(electron_density))
    axes.set_title(time + 'ps')
    axes.plot(x_values, electron_density, label='e$^-$')
    axes.plot(x_values, ion_density, label='ion')
    # axes.set_ylim(1E-6, 1E3)
    # if species == 'ion':
    axes.set_yscale('log')
    axes.set_xlabel('x [$\mu$m]')
    axes.set_ylabel('n$_e$ [cm$^3$]')
    # axes.set_xlim(-50, 5)
    axes.set_ylim(1E14, 1E21)
    # axes.set_xlim(-900, 100)
    axes.set_title(time + 'ns')
    figure.tight_layout()
    plt.legend(frameon=False)
    plt.savefig(folderpath + str(filenumber) + '_centraldensity.png', dpi=200)
    plt.close()
    return x_values, electron_density, ion_density

def central_lineout_mass(x, y, targ, cham, folderpath, filenumber, time, figure, axes):
    if resolution[0] % 2 == 1:
        mid = int((resolution[0] + 1) / 2)
    else:
        mid = int(resolution[0] / 2)

    targ = targ[mid, :]
    cham = cham[mid, :]
    x_values = np.linspace(x[0], x[1], len(targ))
    axes.set_title(time + 'ps')
    axes.plot(x_values, targ, label='target')
    axes.plot(x_values, cham, label='gas')
    # axes.set_ylim(1E-6, 1E3)
    # if species == 'ion':
    axes.set_yscale('log')
    axes.set_xlabel('x [$\mu$m]')
    axes.set_ylabel('mass')
    # axes.set_xlim(-50, 5)
    axes.set_ylim(1E-6, 5)
    # axes.set_xlim(-900, 100)
    axes.set_title(time + 'ps')
    figure.tight_layout()
    plt.legend(frameon=False)
    plt.savefig(folderpath + str(filenumber) + '_centraldensity_mass.png', dpi=200)
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

        x = leftedge[0] - centre[0], rightedge[0] - centre[0]
        y = leftedge[1] - centre[1], rightedge[1] - centre[1]
        
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

    # 2d colorplots
    plot_2D(x, y, electron_temp, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict["tele"])
    plot_2D(x, y, ion_temp, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict["tion"])
    plot_2D(x, y, rad_temp, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict["trad"])
    plot_2D(x, y, electron_density, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict["edens"])
    plot_2D(x, y, ion_density, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict["iondens"])
    if plot_B:
        plot_2D(x, y, Bz, folderpath + folders, hydro_file, time, dictionary=quantity_dict.dict["magz"])

    if counter == len(files)-1:
        make_movie(folderpath + folders, files, 'iondens')
        make_movie(folderpath + folders, files, 'edens')
        make_movie(folderpath + folders, files, 'tele')
        make_movie(folderpath + folders, files, 'tion')
        make_movie(folderpath + folders, files, 'trad')
        if plot_B:
            make_movie(folderpath + folders, files, 'magz')

        files = glob.glob(folderpath + folders + '/*.png')
        for f in files:
            os.remove(f)

    counter +=1




