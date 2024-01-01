import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def get_quantity(frb, quantity):
    data = (frb[quantity].d).T
    data = data[~(data == 0).all(1)]
    return data

def get_arbitrary_line(x, y, data, initial_coordinate, final_coordinate):
    # get the straight line between the 2 points
    slope = (initial_coordinate[1] - final_coordinate[1])/(initial_coordinate[0] - final_coordinate[0])
    intercept = final_coordinate[1] - slope*final_coordinate[0]
    
    x_lineout = np.zeros_like(x)
    y_lineout = np.zeros_like(x)
    data_lineout = np.zeros_like(x)
    # for every x-cordinate, plot the correspodning y cordinate
    x_indices = np.linspace(0, len(x)-1, len(x))
    y_index = np.zeros_like(x)
    for x_index in x_indices:
        # find the corresponding y for every x position according to that equation
        y_plot = slope*x[int(x_index)] + intercept
        corresponding_y_index = np.argmin(np.abs(y-y_plot))
        y_index[int(x_index)] = corresponding_y_index
        x_lineout[int(x_index)] = x[int(x_index)] # this literally does nothing other than it's here for consistency
        y_lineout[int(x_index)] = y[corresponding_y_index]
        data_lineout[int(x_index)] = data[corresponding_y_index, int(x_index)]
    # clip all the data 
    x_index_cut = x_indices[(x_lineout<=initial_coordinate[0]) & (x_lineout>=final_coordinate[0])]
    y_index_cut = y_index[(x_lineout<=initial_coordinate[0]) & (x_lineout>=final_coordinate[0])]
    x_cut = x_lineout[(x_lineout<=initial_coordinate[0]) & (x_lineout>=final_coordinate[0])]
    y_cut = y_lineout[(x_lineout<=initial_coordinate[0]) & (x_lineout>=final_coordinate[0])]
    data_lineout_cut = data_lineout[(x_lineout<=initial_coordinate[0]) & (x_lineout>=final_coordinate[0])]

    return x_index_cut, y_index_cut, data_lineout_cut, x_cut, y_cut

folder = '~/Documents/Hydro/magnetised_shocks/unmag_domenico/unmag28/'
file = 'unmag_shock_1_hdf5_plt_cnt_0020'
filename = folder + file
print(filename)
ds = yt.load(filename)
time = str(round(ds.current_time.d*1E9, 2)) 
slc = ds.slice('z', 0.5)
resolution = [250, 250]

leftedge = ds.domain_left_edge.d*1E4
rightedge = ds.domain_right_edge.d*1E4
centre = [1500, 4900]
target_centre = [(leftedge[0] + rightedge[0]) / 2E4, (leftedge[1] + rightedge[1]) / 2E4, 0] #centre of the simulation in code units (cm)
y = np.linspace(leftedge[0] - centre[0], rightedge[0] - centre[0], resolution[1])
x = np.linspace(leftedge[1] - centre[1], rightedge[1] - centre[1], resolution[0])
frb = slc.to_frb(width=(0.3, 'code_length'), resolution=resolution, height=(0.5, 'code_length'), center=target_centre)

t_ion= get_quantity(frb, 'tion')
t_rad = get_quantity(frb, 'trad')
t_ele = get_quantity(frb, 'tele')
Ye = get_quantity(frb, 'ye')
Sumy = get_quantity(frb, 'sumy')
dens = get_quantity(frb, 'dens')
electron_density = 6.023E23*Ye*dens
ion_density = 6.023E23*Sumy*dens
zbar = Ye/Sumy

# put the initial coordinate in x, y in units of the data (it's been converted to microns at this stage)
# keep the initial coordinate to the right of the final coordinate for now, it's how im clipping the data
initial_coordinate = [-1500, -1000]
final_coordinate = [-3000, 0]

fig, ax = plt.subplots(1, 2, figsize = (10, 6))
cax = ax[0].imshow(electron_density, cmap='Reds', extent=[x[0], x[-1], y[0], y[-1]], norm=colors.LogNorm(vmin=1e14, vmax=1e21), interpolation='nearest', origin='lower', aspect='auto') # log colorbar
# cax = ax.imshow(electron_density, cmap='Reds', extent=[y[0], y[1], x[0], x[-1]], vmin=1, vmax=5000, interpolation='nearest', origin='lower', aspect='auto') # linear colourbar
ax[0].set_xlabel('x [$\mu$m]')
ax[0].set_ylabel('y [$\mu$m]')
cbar = fig.colorbar(cax)
cbar.ax.set_ylabel('n$_e$ [cm$^{-3}$]')
ax[0].set_title(f'{time}ns')
# plt.show()

x_index_cut, y_index_cut, electron_data_lineout, x_cut, y_cut = get_arbitrary_line(x, y, electron_density, initial_coordinate, final_coordinate)
distance = np.sqrt((x_cut-initial_coordinate[0])**2 + (y_cut-initial_coordinate[1])**2)
ax[0].plot(x_cut, y_cut, '--')
ax[1].plot(distance, electron_data_lineout)
ax[1].invert_xaxis() 
ax[1].yaxis.tick_right()
ax[1].yaxis.set_ticks_position('both')
ax[1].set_yscale('log')
ax[1].set_xlabel('distance from initial cooridnate [$\mu$m]')
fig.tight_layout()
plt.show()