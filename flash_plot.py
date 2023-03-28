import numpy as np
import matplotlib.pyplot as plt
import yt
import matplotlib.colors as colors

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
files = range(0, 18)
plot_B = False
plot_bdry = False
plot_mass = True

def get_density_hydro(ds, slc, frb):
    density = (frb['dens'].d).T
    density = density[~(density == 0).all(1)]

    YE = (frb['ye'].d).T
    YE = YE[~(YE == 0).all(1)]

    SUMY = (frb['sumy'].d).T
    SUMY = SUMY[~(SUMY == 0).all(1)]

    Tele = (frb['tele'].d).T
    Tele = Tele[~(Tele == 0).all(1)]

    Tion = (frb['tion'].d).T
    Tion = Tion[~(Tion == 0).all(1)]

    Trad = (frb['trad'].d).T
    Trad = Trad[~(Trad == 0).all(1)]
    return density, YE, SUMY, Tele, Tion, Trad

def get_B_field(ds, slc, frb):
    Bx = (frb['magx'].d).T * np.sqrt(4*np.pi)/1E4
    Bx = Bx[~(Bx == 0).all(1)]

    By = (frb['magy'].d).T* np.sqrt(4*np.pi)/1E4
    By = By[~(By == 0).all(1)]

    Bz = (frb['magz'].d).T* np.sqrt(4*np.pi)/1E4
    Bz = Bz[~(Bz == 0).all(1)]
    return Bx, By, Bz


def get_mass_hydro(ds, slc, frb):
    Targ = (frb['targ'].d).T
    Targ = Targ[~(Targ == 0).all(1)]

    Cham = (frb['cham'].d).T
    Cham = Cham[~(Cham == 0).all(1)]
    return Targ, Cham


def get_bdry(ds, slc, frb):
    bdry = (frb['bdry'].d).T 
    return bdry


def density_plot(x, y, density, folderpath, filenumber, time, species):
    fig, ax = plt.subplots()
    cax = ax.imshow(density, cmap='Blues',
             extent=[x[0], x[1], y[0], y[1]], norm=colors.LogNorm(vmin=1E14, vmax=1E21), interpolation='nearest', origin='lower', aspect='auto')
    # cax = ax.imshow(density, norm=colors.LogNorm(vmin=1E-2, vmax=1000), cmap='Blues_r', interpolation='nearest', origin='lower', aspect='auto')
    ax.set_xlabel('x [$\mu$m]')
    ax.set_title(time + 'ps')
    # ax.set_yticks([-10, -5, 0, 5, 10])
    ax.set_ylabel('y [$\mu$m]')
    cbar = fig.colorbar(cax)
    # ax.set_xlim(-900, 100)
    #ax.set_ylim(-10, 10)
    cbar.ax.set_ylabel('n [cm$^{-3}]$')
    fig.tight_layout()
    plt.savefig(folderpath + filenumber + species + '_density.png', dpi=200)
    plt.close()
    #plt.show()

def B_field_plot(x, y, B_field, folderpath, filenumber, time, vector):
    fig, ax = plt.subplots()
    
    # cax = ax.imshow(B_field, cmap='Blues',
            #  extent=[x[0], x[1], y[0], y[1]], norm=colors.LogNorm(vmin=1E-2, vmax=10), interpolation='nearest', origin='lower', aspect='auto')
    cax = ax.imshow(B_field, cmap='RdBu',
             extent=[x[0], x[1], y[0], y[1]], interpolation='nearest', origin='lower', aspect='auto', vmin=-50, vmax=50)
    ax.set_xlabel('x [$\mu$m]')
    ax.set_title(time + 'ps')
    # ax.set_yticks([-10, -5, 0, 5, 10])
    ax.set_ylabel('y [$\mu$m]')
    cbar = fig.colorbar(cax)
    #ax.set_ylim(-10, 10)
    # ax.set_xlim(-900, 100)
    cbar.ax.set_ylabel('B [T]')
    fig.tight_layout()
    plt.savefig(folderpath + filenumber + vector + '.png', dpi=200)
    plt.close()

def bdry_plot(x, y, bdry, folderpath, filenumber, time, vector):
    fig, ax = plt.subplots()
    
    cax = ax.imshow(bdry, cmap='RdBu',
             extent=[x[0], x[1], y[0], y[1]], interpolation='nearest', origin='lower', aspect='auto', vmin=-1, vmax=1)
    ax.set_xlabel('x [$\mu$m]')
    ax.set_title(time + 'ps')
    #ax.set_ylim(-10, 10)
    # ax.set_yticks([-10, -5, 0, 5, 10])
    ax.set_ylabel('y [$\mu$m]')
    cbar = fig.colorbar(cax)
    # ax.set_xlim(-900, 100)
    cbar.ax.set_ylabel('bdry')
    fig.tight_layout()
    plt.savefig(folderpath + filenumber + vector + '.png', dpi=200)
    plt.close()

def mass_plot(x, y, targ, cham, folderpath, filenumber, time):
    fig, ax = plt.subplots()
    cax = ax.imshow(targ, cmap='Blues',
             extent=[x[0], x[1], y[0], y[1]], norm=colors.LogNorm(vmin=1E-8, vmax=1), interpolation='nearest', origin='lower', aspect='auto')
    ax.set_xlabel('x [$\mu$m]')
    ax.set_title(time + 'ps')
    #ax.set_ylim(-10, 10)
    # ax.set_yticks([-10, -5, 0, 5, 10])
    ax.set_ylabel('y [$\mu$m]')
    cbar = fig.colorbar(cax)
    # ax.set_xlim(-900, 100)
    # cbar.ax.set_ylabel('bdry')
    fig.tight_layout()
    plt.savefig(folderpath + filenumber + 'targ' + '.png', dpi=200)
    plt.close()

    fig, ax = plt.subplots()
    cax = ax.imshow(cham, cmap='Blues',
             extent=[x[0], x[1], y[0], y[1]], norm=colors.LogNorm(vmin=1E-8, vmax=1), interpolation='nearest', origin='lower', aspect='auto')
    ax.set_xlabel('x [$\mu$m]')
    ax.set_title(time + 'ps')
    #ax.set_ylim(-10, 10)
    # ax.set_yticks([-10, -5, 0, 5, 10])
    ax.set_ylabel('y [$\mu$m]')
    cbar = fig.colorbar(cax)
    # ax.set_xlim(-900, 100)
    # cbar.ax.set_ylabel('bdry')
    fig.tight_layout()
    plt.savefig(folderpath + filenumber + 'cham' + '.png', dpi=200)
    plt.close()

def temp_plot(x, y, temp, folderpath, filenumber, time, species):
    fig, ax = plt.subplots()
    cax = ax.imshow(temp, cmap='Reds',norm=colors.LogNorm(vmin=1, vmax=1e5), 
                 extent=[x[0], x[1], y[0], y[1]], interpolation='nearest', origin='lower', aspect='auto')
    #cax = ax.imshow(density, norm=colors.LogNorm(vmin=1E-2, vmax=1000), cmap='Blues_r', interpolation='nearest', origin='lower', aspect='auto')
    ax.set_xlabel('x [$\mu$m]')
    #ax.set_ylim(-10, 10)
    ax.set_title(time + 'ps')
    # ax.set_xlim(-1000, 100)
    # ax.set_yticks([-10, -5, 0, 5, 10])
    ax.set_ylabel('y [$\mu$m]')
    cbar = fig.colorbar(cax)
    cbar.ax.set_ylabel('T [eV]')
    fig.tight_layout()
    plt.savefig(folderpath + filenumber + species + '_temp.png', dpi=200)
    plt.close()
    #plt.show()


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

    density, YE, SUMY, Tele, Tion, Trad = get_density_hydro(ds, slc, frb)

    if plot_B:
        Bx, By, Bz = get_B_field(ds, slc, frb)
        #B_field_plot(x, y, Bx, folderpath + folders, hydro_file, time, vector='Bx')
        #B_field_plot(x, y, By, folderpath + folders, hydro_file, time, vector='By')
        B_field_plot(x, y, Bz, folderpath + folders, hydro_file, time, vector='Bz')
        if counter >0:
            print(np.max(np.max(Bz)))
    if plot_bdry:
        bdry = get_bdry(ds, slc, frb)
        bdry_plot(x, y, bdry, folderpath + folders, hydro_file, time, vector='bdry')

    if plot_mass:
        targ, cham = get_mass_hydro(ds, slc, frb)
        mass_plot(x, y, targ, cham, folderpath + folders, hydro_file, time)
        fig4, ax4 = plt.subplots()
        central_lineout_mass(x, y, targ, cham, folderpath + folders, hydro_file, time, figure=fig4, axes=ax4)
    # eV
    electron_temp = Tele*8.62E-5
    ion_temp = Tion*8.62E-5
    rad_temp = Trad*8.62E-5

    # kelvin
    # electron_temp = Tele
    # ion_temp = Tion
    # rad_temp = Trad

    electron_density = 6.023E23*YE*density
    # ionisation = YE/SUMY
    ion_density = 6.023E23*SUMY*density

    # 2d colorplots
    # density_plot(x, y, YE/SUMY, folderpath + folders, hydro_file, time, species='electron')
    density_plot(x, y, electron_density, folderpath + folders, hydro_file, time, species='electron')
    density_plot(x, y, ion_density, folderpath + folders, hydro_file, time, species='ion')
    temp_plot(x, y, electron_temp, folderpath + folders, hydro_file, time, species='electron')
    temp_plot(x, y, ion_temp, folderpath + folders, hydro_file, time, species='ion')
    temp_plot(x, y, rad_temp, folderpath + folders, hydro_file, time, species='rad')
    # 1d plots
    fig, ax = plt.subplots()
    central_lineout(x, y, electron_density, ion_density, folderpath + folders, hydro_file, time, figure=fig, axes=ax)
    # x_values, ion_lineout = central_lineout(x, y, ion_density, folderpath + folders, hydro_file, species='ion', figure=fig, axes=ax)
    fig3, ax3 = plt.subplots()
    central_lineout_temp(x, y, electron_temp, ion_temp, rad_temp, folderpath + folders, hydro_file, time, figure=fig3, axes=ax3)

    if counter == len(files)-1:
        make_movie(folderpath + folders, files, 'ion_density')
        make_movie(folderpath + folders, files, 'electron_density')
        make_movie(folderpath + folders, files, 'electron_temp')
        make_movie(folderpath + folders, files, 'ion_temp')
        make_movie(folderpath + folders, files, 'rad_temp')
        make_movie(folderpath + folders, files, '_centraltemp')
        make_movie(folderpath + folders, files, '_centraldensity')
        if plot_B:
            #make_movie(folderpath + folders, files, 'Bx')
            #make_movie(folderpath + folders, files, 'By')
            make_movie(folderpath + folders, files, 'Bz')
        if plot_bdry:
            make_movie(folderpath + folders, files, 'bdry')
        if plot_mass:
            make_movie(folderpath + folders, files, 'cham')
            make_movie(folderpath + folders, files, 'targ')
            make_movie(folderpath + folders, files, '_centraldensity_mass')

        files = glob.glob(folderpath + folders + '/*.png')
        for f in files:
            os.remove(f)

    counter +=1




