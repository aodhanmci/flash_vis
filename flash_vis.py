import numpy as np
import matplotlib.pyplot as plt
import yt
import matplotlib.colors as colors
from scipy.ndimage.interpolation import rotate
from scipy import interpolate

from make_hydro_movie import *
from format_hydro_files import *

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.size'] = 14

folderpath = "/Users/aodhan/Documents/HYDRO/Perla/"

filename = "lasslab_hdf5_plt_cnt_"
#cartesian
#centre = [0, 0]
#cylindrical
centre = [70, 50]
resolution = [500, 500]
n_crit = 1.01E21
folders = 'Perla_11/'
files = [10, 20, 40, 60, 80, 100, 120, 140]


def rotateImage(img, angle):
    imgR = rotate(img, angle, reshape=False)
    return imgR


def epoch_inter(hydro_x, hydro_y, density, x_boundaries, y_boundaries, dx, dy):
    n_x = int((x_boundaries[1]-x_boundaries[0])/dx)
    n_y = int((y_boundaries[1]-y_boundaries[0])/dy)

    x_centres = np.linspace(x_boundaries[0], x_boundaries[1], n_x)
    y_centres = np.linspace(y_boundaries[0], y_boundaries[1], n_y)

    hydro_grid_x = np.linspace(hydro_x[0], hydro_x[1], density.shape[0])
    hydro_grid_y = np.linspace(hydro_y[0], hydro_y[1], density.shape[1])

    y_min = np.argmin(np.abs(hydro_grid_y-y_boundaries[0]))
    y_max = np.argmin(np.abs(hydro_grid_y-y_boundaries[1]))

    XI, YI = np.meshgrid(x_centres, y_centres)
    f = interpolate.interp2d(hydro_grid_x, hydro_grid_y[y_min:y_max], density[:, y_min:y_max].T, kind='linear')
    znew = f(x_centres, y_centres)

    return znew, x_centres, y_centres


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
    return density, YE, SUMY, Tele, Tion


def density_plot(x, y, density, folderpath, filenumber, time, species):
    fig, ax = plt.subplots()
    cax = ax.imshow(density, norm=colors.LogNorm(vmin=1E-3, vmax=50), cmap='Blues_r',
                 extent=[x[0], x[1], y[0], y[1]], interpolation='nearest', origin='lower', aspect='auto')
    #cax = ax.imshow(density, norm=colors.LogNorm(vmin=1E-2, vmax=1000), cmap='Blues_r', interpolation='nearest', origin='lower', aspect='auto')
    ax.set_xlabel('x($\mu$m)')
    ax.set_title(time + 'ns')
    # ax.set_yticks([-10, -5, 0, 5, 10])
    ax.set_ylabel('y($\mu$m)')
    cbar = fig.colorbar(cax)
    cbar.ax.set_ylabel('n/n$_c$')
    fig.tight_layout()
    plt.savefig(folderpath + filenumber + species + '_density.png', dpi=200)
    plt.close()
    #plt.show()


def temp_plot(x, y, temp, folderpath, filenumber, species):
    fig, ax = plt.subplots()
    cax = ax.imshow(temp, norm=colors.LogNorm(vmin=1E-2, vmax=1E3), cmap='Blues_r',
                 extent=[x[0], x[1], y[0], y[1]], interpolation='nearest', origin='lower', aspect='auto')
    #cax = ax.imshow(density, norm=colors.LogNorm(vmin=1E-2, vmax=1000), cmap='Blues_r', interpolation='nearest', origin='lower', aspect='auto')
    ax.set_xlabel('x($\mu$m)')
    ax.set_yticks([-10, -5, 0, 5, 10])
    ax.set_ylabel('y($\mu$m)')
    cbar = fig.colorbar(cax)
    cbar.ax.set_ylabel('T$_e$')
    fig.tight_layout()
    plt.savefig(folderpath + filenumber + species + '_temp.png', dpi=200)
    plt.close()
    #plt.show()


def central_lineout_temp(x, y, e_temp, folderpath, filenumber, time, species, figure, axes):
    # if resolution[0] % 2 == 1:
    #     mid = int((resolution[0] + 1) / 2)
    # else:
    #     mid = int(resolution[0] / 2)
    y_values = np.linspace(y[0], y[1], len(e_temp[:, 0]))
    zero_line = np.argmin(np.abs(y_values))
    e_temp = e_temp[zero_line, :]
    x_values = np.linspace(x[0], x[1], len(e_temp))
    axes.plot(x_values, e_temp, label=species)
    # axes.set_ylim(1E-6, 1E3)
    # if species == 'ion':
    axes.set_yscale('log')
    axes.set_xlabel('x($\mu$m)')
    axes.set_ylabel('T$_e$')
    axes.set_xlim(-50, 5)
    axes.set_ylim(1, 100)
    axes.set_title(time + 'ns')
    figure.tight_layout()
    plt.legend(frameon=False)
    plt.savefig(folderpath + str(filenumber) + '_central' + species + 'temp.png', dpi=200)
    plt.close()
    return x_values, e_temp


def central_lineout(x, y, density, folderpath, filenumber, time, species, figure, axes):
    # if resolution[0] % 2 == 1:
    #     mid = int((resolution[0] + 1) / 2)
    # else:
    #     mid = int(resolution[0] / 2)
    line_600nm = np.linspace(1E-2, 10, 100)
    line_600nm_x = 0.6*np.ones(np.shape(line_600nm))
    y_values = np.linspace(y[0], y[1], len(density[:, 0]))
    zero_line = np.argmin(np.abs(y_values))
    density_lineout = density[zero_line, :]
    x_values = np.linspace(x[0], x[1], len(density_lineout))
    axes.plot(x_values, density_lineout, label=species)
    axes.plot(line_600nm_x, line_600nm)
    # axes.set_ylim(1E-6, 1E3)
    # if species == 'ion':
    axes.set_yscale('log')
    axes.set_xlabel('x($\mu$m)')
    axes.set_ylabel('n$_e$/n$_c$')
    axes.set_xlim(-50, 5)
    axes.set_ylim(1E-2, 100)
    axes.set_title(time + 'ns')
    figure.tight_layout()
    plt.legend(frameon=False)
    plt.savefig(folderpath + str(filenumber) + '_central.png', dpi=200)
    plt.close()
    return x_values, density_lineout


def ionisation_lineout(x, y, ionisation, folderpath, filenumber, time, species, figure, axes):
    line_600nm = np.linspace(1E-2, 10, 100)
    line_600nm_x = 0.6*np.ones(np.shape(line_600nm))
    y_values = np.linspace(y[0], y[1], len(ionisation[:, 0]))
    zero_line = np.argmin(np.abs(y_values))
    ionisation_lineout = ionisation[zero_line, :]
    x_values = np.linspace(x[0], x[1], len(ionisation_lineout))
    axes.plot(x_values, ionisation_lineout, label=species)
    axes.plot(line_600nm_x, line_600nm)
    # axes.set_ylim(1E-6, 1E3)
    # if species == 'ion':
    axes.set_yscale('log')
    axes.set_xlabel('x($\mu$m)')
    axes.set_ylabel('n$_e$/n$_c$')
    axes.set_xlim(-50, 5)
    axes.set_ylim(0, 3)
    axes.set_title(time + 'ns')
    figure.tight_layout()
    plt.legend(frameon=False)
    plt.savefig(folderpath + str(filenumber) + '_central_ionisation.png', dpi=200)
    plt.close()
    return x_values, ionisation_lineout


def density_lineout(x, y, density, folderpath, filenumber, time, species, figure, axes):
    line_600nm = np.linspace(1E-2, 10, 100)
    line_600nm_x = 0.6*np.ones(np.shape(line_600nm))
    y_values = np.linspace(y[0], y[1], len(density[:, 0]))
    zero_line = np.argmin(np.abs(y_values))
    mass_density_lineout = density[zero_line, :]
    x_values = np.linspace(x[0], x[1], len(mass_density_lineout))
    axes.plot(x_values, mass_density_lineout, label=species)
    axes.plot(line_600nm_x, line_600nm)
    # axes.set_ylim(1E-6, 1E3)
    # if species == 'ion':
    axes.set_yscale('log')
    axes.set_xlabel('x($\mu$m)')
    axes.set_ylabel('n$_e$/n$_c$')
    axes.set_xlim(-50, 5)
    axes.set_ylim(1E-2, 1)
    axes.set_title(time + 'ns')
    figure.tight_layout()
    plt.legend(frameon=False)
    plt.savefig(folderpath + str(filenumber) + '_density_lineout.png', dpi=200)
    plt.close()
    return x_values, mass_density_lineout

counter=0
for filenumber in files:

    hydro_file = format_file_number(filenumber)
    file = folderpath + folders + filename + hydro_file

    ds = yt.load(file)
    #cartesian
    # slc = ds.slice('z', 0.5)
    #cylindrical
    slc = ds.slice(2, 0)
    time = round(ds.current_time.d*1E9)
    
    time = str(time)
    if counter == 0:
        leftedge = ds.domain_left_edge.d*1E4
        rightedge = ds.domain_right_edge.d*1E4

        x = leftedge[1] - centre[0], rightedge[1] - centre[0]
        y = leftedge[0] - centre[1], rightedge[0] - centre[1]
        
        # cylindrical
        target_centre = [rightedge[0] / 2E4, rightedge[1] / 2E4, 0]
        # cartesian
        # target_centre = [0, 0, 0]
    frb = slc.to_frb(height=(0.01, 'code_length'), resolution=resolution, width=(0.01, 'code_length'),
                         center=target_centre)
    bounds = np.array(frb.bounds)
    x = (bounds[2] * 1E4) - centre[0], (bounds[3]*1E4) - centre[0]
    y = (bounds[0] * 1E4) - centre[1], (bounds[1]*1E4) - centre[1]

    density, YE, SUMY, Tele, Tion = get_density_hydro(ds, slc, frb)
    electron_temp = Tele
    ion_temp = Tion
    # electron_density = 6.023E23*YE*density / n_crit
    ionisation_level = 3.5
    electron_density = 6.023E23*ionisation_level*SUMY*density / n_crit
    
    ionisation = YE/SUMY
    # ion_density = 6.023E23*SUMY*density / n_crit

    #electron_density[:, 515:] = 0
    # 2d colorplots
    # density_plot(x, y, YE/SUMY, folderpath + folders, hydro_file, time, species='electron')
    density_plot(x, y, electron_density, folderpath + folders, hydro_file, time, species='electron')
    #density_plot(x, y, ion_density, folderpath + folders[filenumber], hydro_file, species='ion')
    temp_plot(x, y, electron_temp, folderpath + folders, hydro_file, species='electron')
    temp_plot(x, y, ion_temp, folderpath + folders, hydro_file, species='ion')
    # rotated = rotateImage(electron_density, 33)
    # rot_x = [x[0], x[1]]
    # rot_y = [y[0], y[1]]
    # #rotated[rotated < 1E-3] = 1E-3
    # density_plot(rot_x, rot_y, rotated, folderpath, hydro_file, species='electron_rotate')


    # # 1d lineouts

    # if filenumber == 0:
    #     ion_fig, ion_ax = plt.subplots()
    #     species_name = 'e-_ionised'
    # else:
    #     species_name = 'e-_no_ionisation'

    fig, ax = plt.subplots()
    x_values, electron_lineout = central_lineout(x, y, electron_density, folderpath + folders, hydro_file, time, species='e-', figure=fig, axes=ax)
    fig2, ax2 = plt.subplots()
    x_values, ionisation_line = ionisation_lineout(x, y, ionisation, folderpath + folders, hydro_file, time, species='e-', figure=fig2, axes=ax2)
    # x_values, ion_lineout = central_lineout(x, y, ion_density, folderpath + folders, hydro_file, species='ion', figure=fig, axes=ax)
    fig3, ax3 = plt.subplots()
    x_values, temp_lineout = central_lineout_temp(x, y, electron_temp, folderpath + folders, hydro_file, time, species='electron', figure=fig3, axes=ax3)
    # interpolate onto epoch grid
    fig4, ax4 = plt.subplots()
    x_values, density_line = density_lineout(x, y, density, folderpath + folders, hydro_file, time, species='e-', figure=fig4, axes=ax4)
    
    fig5, ax5 = plt.subplots()
    x_values, ion_temp_lineout = central_lineout_temp(x, y, ion_temp, folderpath + folders, hydro_file, time, species='ion', figure=fig5, axes=ax5)

    # dx = 0.1
    # dy = 0.1
    # x_boundaries = [-70, 5]
    # y_boundaries = [-40, 40]
    #
    # epoch_dens, epoch_x, epoch_y = epoch_inter(x, y, electron_density.T, x_boundaries, y_boundaries, dx, dy)
    # epoch_temp, epoch_x, epoch_y = epoch_inter(x, y, electron_temp.T, x_boundaries, y_boundaries, dx, dy)
    # #  cut the edges off
    # cut_data = new_data[250:1750, :]
    # cut_data[cut_data < 1E-3] = 0
    # plt.figure(3)
    # plt.imshow(epoch_dens, extent=(x_boundaries[0], x_boundaries[1], y_boundaries[0], y_boundaries[1]),
    #             norm=colors.LogNorm(vmin=1E-3, vmax=10), cmap='Blues_r',
    #             interpolation='nearest', origin='lower', aspect='auto')
    # plt.xlim(-70, 5)
    # plt.ylim(-30, 30)
    # plt.show()
    #
    # with open(folderpath + folders + "electron_density.dat", 'wb') as f:
    #     cut_data = epoch_dens
    #     cut_data.tofile(f)


    #     f.close()
    # with open(folderpath + "electron_temp.dat", 'wb') as f:
    #     #epoch_temp[:, 3200:] = 0
    #     epoch_temp[epoch_dens < 1E-3] = 0
    #     for x_ind in range(0, len(epoch_x)):
    #         for y_index in range(0, len(epoch_y)):
    #             if -2 < epoch_x[x_ind] < 2:
    #                 if epoch_temp[y_index, x_ind] < 1:
    #                     epoch_temp[y_index, x_ind] = 1E4
    #     #epoch_temp[any(epoch_dens < 1E-3 and epoch_x < 1E-6)] = 1E4
    #     #epoch_temp[epoch_dens >= 0] = 0
    #     #epoch_temp = epoch_temp[750, :]
    #     epoch_temp.tofile(f)
    #     f.close()
    #     epoch_dens_1D = epoch_dens_1D/1.01E27
    #     plt.plot(epoch_x, epoch_dens_1D)
    #     plt.yscale('log')
    #     plt.xlim(-20, 10)


    # scale_fit_1 = electron_lineout[np.logical_and(x_values > -50, x_values < 0)]
    # energy_fit_1 = x_values[np.logical_and(x_values > -50, x_values < 0)]
    # log_density_1 = np.log(scale_fit_1)

    # scale_fit_2 = electron_lineout[np.logical_and(x_values > -50, x_values < -30)]
    # energy_fit_2 = x_values[np.logical_and(x_values > -50, x_values < -30)]
    # log_density_2 = np.log(scale_fit_2)
    
    # fit_values_1 = np.polyfit(energy_fit_1, log_density_1, 1)
    # fit_values_2 = np.polyfit(energy_fit_2, log_density_2, 1)
    
    # y_fit_1 = np.exp(energy_fit_1*fit_values_1[0]+fit_values_1[1]) 
    # print(fit_values_1[0])
    # y_fit_2 = np.exp(energy_fit_2*fit_values_2[0]+fit_values_2[1]) 
    # plt.plot(x_values, electron_lineout)
    # plt.plot(energy_fit_1, y_fit_1)
    # plt.plot(energy_fit_2, y_fit_2)
    # plt.ylim(1E-3, 5)
    # plt.xlim(-50, 5)
    # plt.yscale('log')
    # plt.show()

    if counter == len(files)-1:
        make_movie(folderpath + folders, files, 'electron_density')
        make_movie(folderpath + folders, files, '_central')
        make_movie(folderpath + folders, files, '_centralelectrontemp')
        make_movie(folderpath + folders, files, '_centraliontemp')
        make_movie(folderpath + folders, files, '_central_ionisation')
        make_movie(folderpath + folders, files, '_density_lineout')
        make_movie(folderpath + folders, files, 'electron_temp')


    #     files = glob.glob(folderpath + '*.png')
    #     #for f in files:
    #         #os.remove(f)

    counter +=1




