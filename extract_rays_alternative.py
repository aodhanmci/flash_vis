import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian(x, a, b, c):
    return a * np.exp(-(x - b)**2 / (c**2))

folder = '/Users/AMcilvenny/Documents/Hydro/magnetised_shocks/unmag_domenico/unmag7/'

# Load the ray data from the HDF5 file
with h5py.File(folder + 'unmag_shock_7_hdf5_plt_cnt_0004', 'r') as f:
    x_ray = np.array(f['RayData'][:, 1])
    y_ray = np.array(f['RayData'][:, 2])
    z_ray = np.array(f['RayData'][:, 3])
    ray_power_ray = np.array(f['RayData'][:, 4]/1e7) # converts from erg/s to joules/s
print(np.nanmin(ray_power_ray))
print(np.nanmax(ray_power_ray))

filter = (ray_power_ray < 3e6) & ((0.15-0.0001 < z_ray) & (z_ray < 0.15+0.0001))
x_ray = x_ray[~filter]
y_ray = y_ray[~filter]
z_ray = z_ray[~filter]
ray_power_ray = ray_power_ray[~filter]
x_lineout = np.linspace(0.1, 0.2, 50)
y_lineout = -np.tan(0*np.pi/180)*x_lineout+0.425

offset = 0.001

# Parallel lines equations
y_upper = y_lineout + offset
y_lower = y_lineout - offset

condition = []
for xi, yi, zi in zip(x_ray, y_ray, z_ray):
    y_on_line = -np.tan(0 * np.pi / 180) * xi + 0.425
    condition.append((y_on_line - offset < yi < y_on_line + offset) and (0.15-0.001 < zi < 0.15+0.001) and (np.min(x_lineout) < xi < np.max(x_lineout))) 

# Filter points
x_filtered, y_filtered, z_filtered, ray_filtered = x_ray[condition], y_ray[condition], z_ray[condition], ray_power_ray[condition]
cross_section = np.sqrt(x_filtered**2 + y_filtered**2)

plt.figure(figsize=(6, 6))
scatter = plt.scatter(x_ray, y_ray, c=ray_power_ray, cmap='inferno', s=3)
plt.plot(x_lineout, y_lineout)
plt.plot(x_lineout, y_upper)
plt.plot(x_lineout, y_lower)
plt.scatter(x_filtered, y_filtered, label="Points within parallel lines")
plt.colorbar(scatter) 
plt.xlabel('X (cm)')
plt.ylabel('Ray Power (W)')
# plt.title('Ray Power for y=0.45')
# plt.xlim(0, 0.3)
# plt.ylim(0, 0.3)
plt.tight_layout()
plt.show()
ray_fwhm_filter = ray_filtered[ray_filtered>2e5]
distance_fwhm_filter = cross_section[ray_filtered>2e5]
plt.figure(figsize=(6, 6))
# plt.scatter(cross_section, ray_filtered)
# Perform the curve fit
initial_guess = [max(ray_fwhm_filter), np.mean(distance_fwhm_filter), np.std(distance_fwhm_filter)]
# Extracting the parameters
popt, pcov = curve_fit(gaussian, distance_fwhm_filter, ray_fwhm_filter, p0=initial_guess)
a, b, c = popt
print(popt)
# Generating y-values based on the fit for plotting
fit_y = gaussian(distance_fwhm_filter, a, b, c)
e_squared = 0.5*a
left_side = distance_fwhm_filter[distance_fwhm_filter < b]
right_side = distance_fwhm_filter[distance_fwhm_filter > b]

# Find the closest points to the 1/e^2 level on each side
left_index = int(np.argmin(np.abs(a * np.exp(-(left_side - b)**2 / (c**2)) - e_squared)))
right_index = int(np.argmin(np.abs(a * np.exp(-(right_side - b)**2 / (c**2)) - e_squared)))
distance = distance_fwhm_filter[int(left_index)], distance_fwhm_filter[int(right_index)]
waist = ((distance[1] - distance[0]))*np.sqrt(2*np.log(2))
print(waist)
plt.scatter(distance_fwhm_filter, ray_fwhm_filter)
plt.scatter(distance_fwhm_filter, fit_y, label='Fitted Gaussian', color='red')
# plt.scatter(distance, [e_squared, e_squared])
plt.show()
print(f'area = {a*2*np.sqrt(2*np.pi)*c}')
intensity = 2*1.35e9/(np.pi*(waist/2)**2)/1e14
print(intensity)