import glob
import os

import numpy as np
from scipy import signal
from scipy.ndimage import zoom
from spectresc import spectres
from matplotlib.pyplot import *

ion()

arc_spec_blue_count = []
arc_spec_red_count = []

arc_spec_blue_coefficients = []
arc_spec_red_coefficients = []

folder_list = next(os.walk(os.path.join("..", "at2021aekl_output")))[1]
for folder in folder_list:
    # science blue
    blue_temp = np.loadtxt(
        os.path.join(
            "..",
            "at2021aekl_output",
            folder,
            "blue_reduced_science_arc_spec.csv",
        ),
        delimiter=",",
    )
    blue_temp -= np.nanpercentile(blue_temp, 2.5)
    arc_spec_blue_count.append(blue_temp)
    # standard blue
    blue_temp = np.loadtxt(
        os.path.join(
            "..",
            "at2021aekl_output",
            folder,
            "blue_reduced_standard_arc_spec.csv",
        ),
        delimiter=",",
    )
    blue_temp -= np.nanpercentile(blue_temp, 2.5)
    arc_spec_blue_count.append(blue_temp)
    # science red
    red_temp = np.loadtxt(
        os.path.join(
            "..",
            "at2021aekl_output",
            folder,
            "red_reduced_science_arc_spec.csv",
        ),
        delimiter=",",
    )
    red_temp -= np.nanpercentile(red_temp, 2.5)
    arc_spec_red_count.append(red_temp)
    # standard red
    red_temp = np.loadtxt(
        os.path.join(
            "..",
            "at2021aekl_output",
            folder,
            "red_reduced_standard_arc_spec.csv",
        ),
        delimiter=",",
    )
    red_temp -= np.nanpercentile(red_temp, 2.5)
    arc_spec_red_count.append(red_temp)


arc_spec_blue_shift = []
shifted_blue_arcspec = []
for i, arcspec in enumerate(arc_spec_blue_count[1:]):
    con_min = len(arcspec) // 2 - 500
    con_max = len(arcspec) // 2 + 500
    convolved = signal.correlate(
        zoom(arc_spec_blue_count[0], 10.0), zoom(arcspec, 10.0), mode="same"
    )
    shift_blue = (
        np.argmax(convolved[con_min * 10 : con_max * 10]) - 5000
    ) / 10.0
    arc_spec_blue_shift.append(shift_blue)
    shifted_blue_arcspec.append(
        spectres(
            np.arange(len(arcspec)) - shift_blue,
            np.arange(len(arcspec)),
            arcspec,
        )
    )


arc_spec_red_shift = []
shifted_red_arcspec = []
for i, arcspec in enumerate(arc_spec_red_count[1:]):
    con_min = len(arcspec) // 2 - 500
    con_max = len(arcspec) // 2 + 500
    convolved = signal.correlate(
        zoom(arc_spec_red_count[0], 10.0), zoom(arcspec, 10.0), mode="same"
    )
    shift_red = (
        np.argmax(convolved[con_min * 10 : con_max * 10]) - 5000
    ) / 10.0
    arc_spec_red_shift.append(shift_red)
    shifted_red_arcspec.append(
        spectres(
            np.arange(len(arcspec)) - shift_red,
            np.arange(len(arcspec)),
            arcspec,
        )
    )


figure(2)
clf()
plot(np.array(shifted_blue_arcspec).T)

figure(3)
clf()
plot(np.array(shifted_red_arcspec).T)
