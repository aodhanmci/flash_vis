import pathlib
import cv2
from cv2 import VideoWriter, VideoWriter_fourcc
import numpy as np
from format_hydro_files import *
import os
import glob


def make_movie(folderpath, files, quantity):
    initialmoviefile = format_file_number(files[0])
    initial_file_path = folderpath + str(initialmoviefile) + quantity + '.png'
    print(initial_file_path)
    h = cv2.imread(initial_file_path)
    height = h.shape[0]
    width = h.shape[1]
    FPS = 1.0
    fourcc = VideoWriter_fourcc(*'XVID')
    video = VideoWriter(folderpath + quantity + '.mp4', fourcc, float(FPS), (width, height))
    for filenumber in files:
        hydro_file = format_file_number(filenumber)
        filename = folderpath + hydro_file + quantity + '.png'
        filepath = pathlib.Path(filename)
        if filepath.exists():
            h = cv2.imread(filename)
            video.write(np.uint8(h))
        else:
             print(filename + ' does not exist')
    video.release()
