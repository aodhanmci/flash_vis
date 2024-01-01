
dict_2D = {
    "tele":{"name":'tele', "vmin":1, "vmax":10000, "cmap":'Reds', "cbar_label":'T$_{e}$ [eV]'},
    "tion":{"name":'tion', "vmin":1, "vmax":10000, "cmap":'Reds', "cbar_label":'T$_{i}$ [eV]'},
    "trad":{"name":'trad', "vmin":1, "vmax":500, "cmap":'Reds', "cbar_label":'T$_{rad}$ [eV]'},
    "edens_cm3":{"name":'edens', "vmin":1e13, "vmax":1e18, "cmap":'Blues', "cbar_label":'n$_e$ [cm$^{-3}$]'},
    "iondens_cm3":{"name":'iondens', "vmin":1e13, "vmax":1e18, "cmap":'Blues', "cbar_label":'n$_i$ [cm$^{-3}$]'},
    "edens":{"name":'edens', "vmin":1e-3, "vmax":10, "cmap":'Blues', "cbar_label":'n$_e$ [n$_c$]'},
    "iondens":{"name":'iondens', "vmin":1e-3, "vmax":10, "cmap":'Blues', "cbar_label":'n$_i$ [n$_c$]'},
    "magz":{"name":'magz', "vmin":0, "vmax":10, 'cmap':'RdBu', "cbar_label":'B$_z$ [T]'}}

dict_1D = {
    "density_cm3":{"min":1E12, "max":1E21, "ylabel":'n [cm$^{-3}$]'},
    "density":{"min":1e-2, "max":10, "ylabel":'n [n$_c$]'},
    "temp_eV":{"min":1, "max":10000, "ylabel":'T [eV]'},
    "temp_K":{"min":1e5, "max":1e9, "ylabel":'T [K]'}
}