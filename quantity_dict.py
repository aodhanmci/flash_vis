
dict_2D = {
    "tele":{"name":'tele', "vmin":1, "vmax":5e5, "cmap":'Reds', "cbar_label":'T$_{e}$ [eV]'},
    "tion":{"name":'tion', "vmin":1, "vmax":5e5, "cmap":'Reds', "cbar_label":'T$_{i}$ [eV]'},
    "trad":{"name":'trad', "vmin":1, "vmax":5e5, "cmap":'Reds', "cbar_label":'T$_{rad}$ [eV]'},
    "edens":{"name":'edens', "vmin":1e14, "vmax":1e21, "cmap":'Blues', "cbar_label":'n$_e$ [cm$^{-3}$]'},
    "iondens":{"name":'iondens', "vmin":1e14, "vmax":1e21, "cmap":'Blues', "cbar_label":'n$_i$ [cm$^{-3}$]'},
    "magz":{"name":'magz', "vmin":0, "vmax":10, 'cmap':'RdBu', "cbar_label":'B$_z$ [T]'}}

dict_1D = {
    "density":{"min":1E14, "max":1E21, "ylabel":'n [cm$^{-3}$]'},
    "temp_eV":{"min":1, "max":1e4, "ylabel":'T [eV]'},
    "temp_K":{"min":1e5, "max":1e9, "ylabel":'T [K]'}
}