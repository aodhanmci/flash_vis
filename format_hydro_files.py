
def format_file_number(filenumber):
    if filenumber < 10:
        hydro_file = '000' + str(filenumber)
    elif filenumber < 100:
        hydro_file = '00' + str(filenumber)
    elif filenumber < 1000:
        hydro_file = '0' + str(filenumber)
    else:
        hydro_file = str(filenumber)
    return hydro_file
