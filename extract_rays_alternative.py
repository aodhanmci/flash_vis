import h5py
import numpy as np
import matplotlib.pyplot as plt

# Load the ray data from the HDF5 file
with h5py.File('unmag_shock_1_hdf5_plt_cnt_0015', 'r') as f:
    x_ray = f['RayData'][:, 1]
    y_ray = f['RayData'][:, 2]
    ray_power_ray = f['RayData'][:, -1]

# Extract the ray power for rays with y values close to 0.45
mask_y_045 = np.isclose(y_ray, 0.45, atol=1e-4)
x_045 = x_ray[mask_y_045]
ray_power_045 = ray_power_ray[mask_y_045]

mask_x_filter = ray_power_045 > 1.6e16
x_045_filtered = x_045[mask_x_filter]
ray_power_045_filtered = ray_power_045[mask_x_filter]
ray_FWHM = ray_power_045_filtered[ray_power_045_filtered > np.max(ray_power_045_filtered)/2]
x_FWHM = x_045_filtered[ray_power_045_filtered > np.max(ray_power_045_filtered)/2]
FWHM = (np.max(x_FWHM)-np.min(x_FWHM))/np.sqrt(2)
print(f'FWHM:{FWHM*1e4}')
waist = FWHM/1.177
print(f'waist:{waist*1e4}')
print(np.max(ray_power_045_filtered)/(2*waist))
plt.figure(figsize=(10, 6))
plt.plot(x_045_filtered, ray_power_045_filtered, '.')
plt.xlabel('X (cm)')
plt.ylabel('Ray Power (W)')
plt.title('Ray Power for y=0.45')
plt.tight_layout()
plt.grid(True)
plt.show()
