import h5py

def list_datasets(hdf5_group, prefix=''):
    """
    Recursively lists all datasets in an HDF5 group.
    :param hdf5_group: Group object to start searching from
    :param prefix: String to prefix to group/dataset names (used in recursion)
    """
    for key in hdf5_group.keys():
        item = hdf5_group[key]
        path = f"{prefix}/{key}"
        
        if isinstance(item, h5py.Dataset):
            print(f"Dataset: {path}")
        elif isinstance(item, h5py.Group):
            print(f"Group: {path}")
            list_datasets(item, path)
folder = '/Users/AMcilvenny/Documents/Hydro/magnetised_shocks/unmag_domenico/unmag7/unmag_shock_7_hdf5_plt_cnt_0002'
# Replace 'your_file.h5' with the path to your HDF5 file
with h5py.File(folder, 'r') as file:
    list_datasets(file)