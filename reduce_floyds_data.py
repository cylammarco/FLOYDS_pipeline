#!/usr/bin/env python3

import copy
import datetime
import logging
import os
import sys

from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
from astroscrappy import detect_cosmics
from aspired import image_reduction
from aspired import spectral_reduction
from ccdproc import Combiner
from matplotlib import pyplot as plt
from scipy.signal import medfilt
from spectres import spectres
from statsmodels.nonparametric.smoothers_lowess import lowess
import numpy as np
import yaml


# Line list
atlas_Hg_red = [5460.7348, 5769.5982, 5790.6630]
atlas_Ar_red = [
    6965.4307,
    7067.2175,
    7147.0416,
    7272.9359,
    7383.9805,
    7503.8691,
    7635.1056,
    7723.7599,
    7948.1764,
    8014.7857,
    8115.3108,
    8264.5225,
    8424.6475,
    8521.4422,
    9122.9674,
    9224.4992,
    9657.7863,
]
element_Hg_red = ["Hg"] * len(atlas_Hg_red)
element_Ar_red = ["Ar"] * len(atlas_Ar_red)

atlas_Ar_blue = [4158.590, 4200.674, 4272.169, 4300.101, 4510.733]
atlas_Hg_blue = [
    3650.153,
    4046.563,
    4077.8314,
    4358.328,
    5460.7348,
    5769.598,
    5790.663,
]
atlas_Zn_blue = [4722.1569]
element_Hg_blue = ["Hg"] * len(atlas_Hg_blue)
element_Ar_blue = ["Ar"] * len(atlas_Ar_blue)
element_Zn_blue = ["Zn"] * len(atlas_Zn_blue)

# Set the frame
red_spatial_mask = np.arange(0, 330)
blue_spatial_mask = np.arange(335, 512)
red_spec_mask = np.arange(0, 1800)
blue_spec_mask = np.arange(500, 2050)


def flux_diff(ratio, a, b):
    diff = a * ratio - b
    mask = diff < np.nanpercentile(diff, 95)
    return np.nansum(diff[mask] ** 2.0)


# Apply frindg correction onto a target_onedspec
def fringe_correction(target, flat):
    flat_count = flat.spectrum_list[0].count
    flat_continuum = lowess(
        flat_count, np.arange(len(flat_count)), frac=0.04, return_sorted=False
    )
    flat_continuum_divided = flat_count / flat_continuum
    flat_continuum_divided[:500] = 1.0
    # Apply the flat correction
    target.spectrum_list[0].count /= flat_continuum_divided


HERE = os.getcwd()

# Get config file
try:
    param_filename = sys.argv[1]
    params_path = os.path.join(HERE, param_filename)
except Exception as e:
    print(e)
    params_path = os.path.join(HERE, "floyds_2022jdf_20220629_2932516.yaml")

if not os.path.isabs(params_path):
    params_path = os.path.abspath(params_path)

print("Reading parameters from " + params_path + ".")

with open(params_path, "r") as stream:
    params = yaml.safe_load(stream)

# Set up the logger
logger_name = "floyds_reduction_" + params["target_name"]
logger = logging.getLogger(logger_name)
log_level = params["log_level"].lower()
if log_level == "critical":
    logger.setLevel(logging.CRITICAL)
elif log_level == "error":
    logger.setLevel(logging.ERROR)
elif log_level == "warning":
    logger.setLevel(logging.WARNING)
elif log_level == "info":
    logger.setLevel(logging.INFO)
elif log_level == "debug":
    logger.setLevel(logging.DEBUG)

# log file output filename
log_file_name = params["log_file_name"]
# if log file name is set to 'default', use the logger name and current time
if log_file_name is None:
    log_file_name = "{}_{}.log".format(
        logger_name,
        datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S"),
    )

# log file output folder
log_file_folder = params["log_file_folder"]
params_folder_path = os.path.dirname(params_path)
output_folder_path = os.path.join(
    params_folder_path, params["output_folder_path"]
)

if not os.path.exists(output_folder_path):
    os.makedirs(output_folder_path)

if log_file_folder is None:
    log_file_folder = output_folder_path

# configure the logger display format
formatter = logging.Formatter(
    "[%(asctime)s] %(levelname)s [%(filename)s:%(lineno)d] " "%(message)s",
    datefmt="%a, %d %b %Y %H:%M:%S",
)

handler = logging.FileHandler(
    os.path.join(log_file_folder, log_file_name), "a+"
)
handler.setFormatter(formatter)
logger.addHandler(handler)

# Get the working and output directory
logger.info("The parameter file is at: {}".format(params_folder_path))

# If the output folder path is default use the logger name.
if output_folder_path is None:
    output_folder_path = os.path.join(HERE, os.extsep(log_file_name)[0])

if not os.path.exists(output_folder_path):
    os.mkdir(output_folder_path)
    logger.info(
        "The output folder does not exist, a folder is created at: "
        "{}".format(output_folder_path)
    )

logger.info("The output folder is at: {}".format(output_folder_path))

cosmicray = params["cosmicray"]
sigclip = params["sigclip"]
fsmode = params["fsmode"]
psfmodel = params["psfmodel"]
psfsize = params["psfsize"]

if params["hemisphere"] == "north":

    rec_coeff_blue_science = params["science_blue_rectification_coeff_north"]
    rec_coeff_red_science = params["science_red_rectification_coeff_north"]
    rec_coeff_blue_standard = params["standard_blue_rectification_coeff_north"]
    rec_coeff_red_standard = params["standard_red_rectification_coeff_north"]

elif params["hemisphere"] == "south":

    rec_coeff_blue_science = params["science_blue_rectification_coeff_south"]
    rec_coeff_red_science = params["science_red_rectification_coeff_south"]
    rec_coeff_blue_standard = params["standard_blue_rectification_coeff_south"]
    rec_coeff_red_standard = params["standard_red_rectification_coeff_south"]

else:

    print(
        "Only north or south are acceptable. You have provided {}.".format(
            params["hemisphere"]
        )
    )
#
#
#
# Start data reduction here
#
#
#
extracted_twodspec = {
    "science": {
        "red": spectral_reduction.TwoDSpec(
            log_file_name=os.path.join(log_file_folder, log_file_name),
        ),
        "blue": spectral_reduction.TwoDSpec(
            log_file_name=os.path.join(log_file_folder, log_file_name),
        ),
    },
    "standard": {
        "red": spectral_reduction.TwoDSpec(
            log_file_name=os.path.join(log_file_folder, log_file_name),
        ),
        "blue": spectral_reduction.TwoDSpec(
            log_file_name=os.path.join(log_file_folder, log_file_name),
        ),
    },
}
img = {
    "science": None,
    "standard": None,
}
ap_trace_position = {}
ap_trace_sigma = {}

target_name = {}
total_exposure_time = 0.0

# Do standard first, so that the ap-trace can be used for curvature correction
# in both science and standard
for frame_type in ["standard", "science"]:

    logger.info("Start working on the {} frames.".format(frame_type))

    # if the path to the light frame is provided
    if params[frame_type + "_light_frame"] is not None:

        light_path = [
            os.path.join(HERE, i) for i in params[frame_type + "_light_frame"]
        ]

    # If the paths to the frames are not provided, search into the folder
    else:

        # Get all the parameters for the file locations
        light_folder = os.path.join(HERE, params[frame_type + "_light_folder"])
        logger.info("Looking for light frames at {}.".format(light_folder))

        light_path = os.listdir(light_folder)
        logger.info("These are the light frames found: {}".format(light_path))

        # If light or flat folders are empty, fail right away
        if light_path == []:

            logger.error(light_folder + " cannot be empty.")

    # if the path to the flat frame is provided
    if params[frame_type + "_light_frame"] is not None:

        flat_path = [
            os.path.join(HERE, i) for i in params[frame_type + "_flat_frame"]
        ]

    else:

        flat_folder = os.path.join(HERE, params[frame_type + "_flat_folder"])
        logger.info("Looking for flat frames at {}.".format(flat_folder))

        flat_path = os.listdir(flat_folder)
        logger.info("These are the flat frames found: {}".format(flat_path))

        if flat_path == []:

            logger.error(flat_folder + " cannot be empty.")

    # if the path to the arc frame is provided
    if params[frame_type + "_light_frame"] is not None:

        arc_path = [
            os.path.join(HERE, i) for i in params[frame_type + "_arc_frame"]
        ]

    else:

        arc_folder = os.path.join(HERE, params[frame_type + "_arc_folder"])
        logger.info("Looking for arc frames at {}.".format(arc_folder))

        arc_path = os.listdir(arc_folder)
        logger.info("These are the arc frames found: {}".format(arc_path))

        # Throw warning with arcs cannot be found.
        if arc_path == []:

            logger.error(arc_folder + " cannot be empty.")

    light_extension = params[frame_type + "_light_extension"]
    logger.info(
        "Looking for light frames with extension {}.".format(light_extension)
    )

    flat_extension = params[frame_type + "_flat_extension"]
    logger.info(
        "Looking for flat frames with extension {}.".format(flat_extension)
    )

    arc_extension = params[frame_type + "_arc_extension"]
    logger.info(
        "Looking for arc frames with extension {}.".format(arc_extension)
    )

    # Get only the light files with the right extension
    light_sigma_clipping = params[frame_type + "_light_sigma_clipping"]
    light_clip_low = params[frame_type + "_light_clip_low"]
    light_clip_high = params[frame_type + "_light_clip_high"]
    light_combine_type = params[frame_type + "_light_combine_type"]
    light_path_new = []
    for i, path_i in enumerate(light_path):
        if not os.path.isabs(path_i):
            path_i = os.path.abspath(os.path.join(light_folder, path_i))
        if os.path.splitext(path_i.lower())[-1][1:] in light_extension:
            if os.path.exists(path_i):
                light_path_new.append(path_i)
            else:
                raise ValueError(path_i + " does not exist.")

    light_path = light_path_new
    logger.info(
        "These are the light frames found, in full path: {}".format(light_path)
    )

    # Get only the flat files with the right extension
    flat_sigma_clipping = params[frame_type + "_flat_sigma_clipping"]
    flat_clip_low = params[frame_type + "_flat_clip_low"]
    flat_clip_high = params[frame_type + "_flat_clip_high"]
    flat_combine_type = params[frame_type + "_flat_combine_type"]
    flat_path_new = []
    for i, path_i in enumerate(flat_path):
        if not os.path.isabs(path_i):
            path_i = os.path.abspath(os.path.join(flat_folder, path_i))
        if os.path.splitext(path_i.lower())[-1][1:] in flat_extension:
            if os.path.exists(path_i):
                flat_path_new.append(path_i)
            else:
                raise ValueError(path_i + " does not exist.")

    flat_path = flat_path_new
    logger.info(
        "These are the flat frames found, in full path: {}".format(flat_path)
    )

    # Get only the arc files with the right extension
    arc_sigma_clipping = params[frame_type + "_arc_sigma_clipping"]
    arc_clip_low = params[frame_type + "_arc_clip_low"]
    arc_clip_high = params[frame_type + "_arc_clip_high"]
    arc_combine_type = params[frame_type + "_arc_combine_type"]
    arc_path_new = []
    for i, path_i in enumerate(arc_path):
        if not os.path.isabs(path_i):
            path_i = os.path.abspath(os.path.join(arc_folder, path_i))
        if os.path.splitext(path_i.lower())[-1][1:] in arc_extension:
            if os.path.exists(path_i):
                arc_path_new.append(path_i)
            else:
                raise ValueError(path_i + " does not exist.")

    arc_path = arc_path_new
    logger.info(
        "These are the arc frames found, in full path: {}".format(arc_path)
    )

    img[frame_type] = image_reduction.ImageReduction(
        log_file_name=os.path.join(log_file_folder, log_file_name),
    )
    logger.info("ImageReduction object created.")

    # Add the light frames
    for i, path_i in enumerate(light_path):
        try:
            light_temp = fits.open(path_i)["SCI"]
        except:
            light_temp = fits.open(path_i)[1]
        light_temp_data = np.array(light_temp.data).astype("float") / float(
            light_temp.header["GAIN"]
        )
        light_temp_cleaned = detect_cosmics(
            light_temp_data,
            gain=float(light_temp.header["GAIN"]),
            readnoise=float(light_temp.header["RDNOISE"]),
            sigclip=sigclip,
            fsmode=fsmode,
            psfmodel=psfmodel,
            psfsize=psfsize,
        )[1]
        img[frame_type].add_light(
            light_temp_cleaned, light_temp.header, light_temp.header["EXPTIME"]
        )
        logger.info("light frame {} added.".format(path_i))
        if frame_type == "science":
            total_exposure_time += float(light_temp.header["EXPTIME"])

    # Add the arc frames
    for i, path_i in enumerate(arc_path):
        try:
            arc_temp = fits.open(path_i)["SCI"]
        except:
            arc_temp = fits.open(path_i)[1]
        img[frame_type].add_arc(arc_temp.data, arc_temp.header)
        logger.info("arc frame {} added.".format(path_i))

    # There is nothing to reduce, but this put things in the right place
    img[frame_type].reduce()
    logger.info("{} image reduced.".format(frame_type))

    # Collect the flat frames
    flat_header = []
    flat_data = []
    for i, path_i in enumerate(flat_path):
        try:
            flat_temp = fits.open(path_i)["SCI"]
        except:
            flat_temp = fits.open(path_i)[1]
        flat_header.append(flat_temp.header)
        flat_data.append(
            CCDData(np.array(flat_temp.data).astype("float"), unit=u.ct)
        )
        logger.info("flat frame {} added.".format(path_i))

    flat_combiner = Combiner(flat_data)
    flat_combiner.sigma_clipping(
        low_thresh=flat_clip_low, high_thresh=flat_clip_high, func=np.ma.median
    )

    if flat_combine_type == "median":

        flat_combined_CCDdata = flat_combiner.median_combine()

    if flat_combine_type in ["average", "mean"]:

        flat_combined_CCDdata = flat_combiner.average_combine()

    # Red & Blue arm
    for arm in ["red", "blue"]:

        if arm == "red":

            spatial_mask = red_spatial_mask
            spec_mask = red_spec_mask
            if frame_type == "science":
                coeff = rec_coeff_red_science
            if frame_type == "standard":
                coeff = rec_coeff_red_standard

        else:

            spatial_mask = blue_spatial_mask
            spec_mask = blue_spec_mask
            if frame_type == "science":
                coeff = rec_coeff_blue_science
            if frame_type == "standard":
                coeff = rec_coeff_blue_standard

        twodspec = extracted_twodspec[frame_type][arm]

        twodspec.add_data(img[frame_type])
        twodspec.set_properties(
            readnoise=float(light_temp.header["RDNOISE"]),
            gain=float(light_temp.header["GAIN"]),
            spatial_mask=spatial_mask,
            spec_mask=spec_mask,
        )
        twodspec.set_readnoise()
        twodspec.set_gain()
        twodspec.set_seeing()
        twodspec.set_exptime()
        twodspec.set_airmass()

        # Trim the arc the same way
        twodspec.apply_mask_to_arc()

        # Get the trace to rectify the image
        if frame_type == "standard":

            twodspec.ap_trace(
                nspec=params[frame_type + "_" + arm + "_aptrace_nspec"],
                smooth=params[frame_type + "_" + arm + "_aptrace_smooth"],
                nwindow=params[frame_type + "_" + arm + "_aptrace_nwindow"],
                spec_sep=params[frame_type + "_" + arm + "_aptrace_spec_sep"],
                trace_width=params[
                    frame_type + "_" + arm + "_aptrace_trace_width"
                ],
                resample_factor=params[
                    frame_type + "_" + arm + "_aptrace_resample_factor"
                ],
                rescale=params[frame_type + "_" + arm + "_aptrace_rescale"],
                scaling_min=params[
                    frame_type + "_" + arm + "_aptrace_scaling_min"
                ],
                scaling_max=params[
                    frame_type + "_" + arm + "_aptrace_scaling_max"
                ],
                scaling_step=params[
                    frame_type + "_" + arm + "_aptrace_scaling_step"
                ],
                percentile=params[
                    frame_type + "_" + arm + "_aptrace_percentile"
                ],
                shift_tol=params[
                    frame_type + "_" + arm + "_aptrace_shift_tol"
                ],
                fit_deg=params[frame_type + "_" + arm + "_aptrace_fit_deg"],
                ap_faint=params[frame_type + "_" + arm + "_aptrace_ap_faint"],
                display=params[frame_type + "_" + arm + "_aptrace_display"],
                renderer=params[frame_type + "_" + arm + "_aptrace_renderer"],
                width=params[frame_type + "_" + arm + "_aptrace_width"],
                height=params[frame_type + "_" + arm + "_aptrace_height"],
                return_jsonstring=params[
                    frame_type + "_" + arm + "_aptrace_return_jsonstring"
                ],
                save_fig=params[frame_type + "_" + arm + "_aptrace_save_fig"],
                fig_type=params[frame_type + "_" + arm + "_aptrace_fig_type"],
                filename=os.path.join(
                    output_folder_path,
                    params[frame_type + "_" + arm + "_aptrace_filename"]
                    + "_"
                    + frame_type,
                ),
                open_iframe=params[
                    frame_type + "_" + arm + "_aptrace_open_iframe"
                ],
            )

            ap_trace_position[arm] = copy.deepcopy(
                twodspec.spectrum_list[0].trace
            )
            ap_trace_sigma[arm] = copy.deepcopy(
                twodspec.spectrum_list[0].trace_sigma
            )

        if frame_type == "science":

            twodspec.add_trace(
                ap_trace_position[arm], ap_trace_sigma[arm], spec_id=0
            )

        # Get the rectification polynomial
        twodspec.get_rectification(
            upsample_factor=params[
                frame_type + "_" + arm + "_rectification_upsample_factor"
            ],
            bin_size=params[
                frame_type + "_" + arm + "_rectification_bin_size"
            ],
            n_bin=params[frame_type + "_" + arm + "_rectification_n_bin"],
            spline_order=params[
                frame_type + "_" + arm + "_rectification_spline_order"
            ],
            order=params[frame_type + "_" + arm + "_rectification_order"],
            coeff=coeff,
            use_arc=params[frame_type + "_" + arm + "_rectification_use_arc"],
            apply=params[frame_type + "_" + arm + "_rectification_apply"],
            display=params[frame_type + "_" + arm + "_rectification_display"],
            renderer=params[
                frame_type + "_" + arm + "_rectification_renderer"
            ],
            width=params[frame_type + "_" + arm + "_rectification_width"],
            height=params[frame_type + "_" + arm + "_rectification_height"],
            return_jsonstring=params[
                frame_type + "_" + arm + "_rectification_return_jsonstring"
            ],
            save_fig=params[
                frame_type + "_" + arm + "_rectification_save_fig"
            ],
            fig_type=params[
                frame_type + "_" + arm + "_rectification_fig_type"
            ],
            filename=os.path.join(
                output_folder_path,
                params[frame_type + "_" + arm + "_rectification_filename"]
                + "_"
                + frame_type,
            ),
            open_iframe=params[
                frame_type + "_" + arm + "_rectification_open_iframe"
            ],
        )

        # Apply the rectification
        twodspec.apply_rectification()

        # Need to store the traces for fringe correction before overwriting them
        # with the new traces
        trace = copy.deepcopy(twodspec.spectrum_list[0].trace)
        trace_sigma = copy.deepcopy(twodspec.spectrum_list[0].trace_sigma)

        # Get the trace again for the rectified image and then extract
        twodspec.ap_trace(
            nspec=params[frame_type + "_" + arm + "_aptrace_nspec"],
            smooth=params[frame_type + "_" + arm + "_aptrace_smooth"],
            nwindow=params[frame_type + "_" + arm + "_aptrace_nwindow"],
            spec_sep=params[frame_type + "_" + arm + "_aptrace_spec_sep"],
            trace_width=params[
                frame_type + "_" + arm + "_aptrace_trace_width"
            ],
            resample_factor=params[
                frame_type + "_" + arm + "_aptrace_resample_factor"
            ],
            rescale=params[frame_type + "_" + arm + "_aptrace_rescale"],
            scaling_min=params[
                frame_type + "_" + arm + "_aptrace_scaling_min"
            ],
            scaling_max=params[
                frame_type + "_" + arm + "_aptrace_scaling_max"
            ],
            scaling_step=params[
                frame_type + "_" + arm + "_aptrace_scaling_step"
            ],
            percentile=params[frame_type + "_" + arm + "_aptrace_percentile"],
            shift_tol=params[frame_type + "_" + arm + "_aptrace_shift_tol"],
            fit_deg=2,
            ap_faint=params[frame_type + "_" + arm + "_aptrace_ap_faint"],
            display=params[frame_type + "_" + arm + "_aptrace_display"],
            renderer=params[frame_type + "_" + arm + "_aptrace_renderer"],
            width=params[frame_type + "_" + arm + "_aptrace_width"],
            height=params[frame_type + "_" + arm + "_aptrace_height"],
            return_jsonstring=params[
                frame_type + "_" + arm + "_aptrace_return_jsonstring"
            ],
            save_fig=params[frame_type + "_" + arm + "_aptrace_save_fig"],
            fig_type=params[frame_type + "_" + arm + "_aptrace_fig_type"],
            filename=os.path.join(
                output_folder_path,
                params[frame_type + "_" + arm + "_aptrace_filename"]
                + "_"
                + frame_type,
            ),
            open_iframe=params[
                frame_type + "_" + arm + "_aptrace_open_iframe"
            ],
        )

        twodspec.ap_extract(
            apwidth=params[frame_type + "_" + arm + "_extract_apwidth"],
            skysep=params[frame_type + "_" + arm + "_extract_skysep"],
            skywidth=params[frame_type + "_" + arm + "_extract_skywidth"],
            skydeg=params[frame_type + "_" + arm + "_extract_skydeg"],
            sky_sigma=params[frame_type + "_" + arm + "_extract_sky_sigma"],
            optimal=params[frame_type + "_" + arm + "_extract_optimal"],
            algorithm=params[frame_type + "_" + arm + "_extract_algorithm"],
            model=params[frame_type + "_" + arm + "_extract_model"],
            lowess_frac=params[
                frame_type + "_" + arm + "_extract_lowess_frac"
            ],
            lowess_it=params[frame_type + "_" + arm + "_extract_lowess_it"],
            lowess_delta=params[
                frame_type + "_" + arm + "_extract_lowess_delta"
            ],
            tolerance=params[frame_type + "_" + arm + "_extract_tolerance"],
            cosmicray_sigma=params[
                frame_type + "_" + arm + "_extract_cosmicray_sigma"
            ],
            max_iter=params[frame_type + "_" + arm + "_extract_max_iter"],
            forced=params[frame_type + "_" + arm + "_extract_forced"],
            variances=params[frame_type + "_" + arm + "_extract_variances"],
            npoly=params[frame_type + "_" + arm + "_extract_npoly"],
            polyspacing=params[
                frame_type + "_" + arm + "_extract_polyspacing"
            ],
            pord=params[frame_type + "_" + arm + "_extract_pord"],
            qmode=params[frame_type + "_" + arm + "_extract_qmode"],
            nreject=params[frame_type + "_" + arm + "_extract_nreject"],
            display=params[frame_type + "_" + arm + "_extract_display"],
            renderer=params[frame_type + "_" + arm + "_extract_renderer"],
            width=params[frame_type + "_" + arm + "_extract_width"],
            height=params[frame_type + "_" + arm + "_extract_height"],
            return_jsonstring=params[
                frame_type + "_" + arm + "_extract_return_jsonstring"
            ],
            save_fig=params[frame_type + "_" + arm + "_extract_save_fig"],
            fig_type=params[frame_type + "_" + arm + "_extract_fig_type"],
            filename=os.path.join(
                output_folder_path,
                params[frame_type + "_" + arm + "_extract_filename"]
                + "_"
                + frame_type,
            ),
            open_iframe=params[
                frame_type + "_" + arm + "_extract_open_iframe"
            ],
        )

        twodspec.extract_arc_spec(
            spec_width=params[frame_type + "_" + arm + "_arc_spec_spec_width"],
            display=params[frame_type + "_" + arm + "_arc_spec_display"],
            renderer=params[frame_type + "_" + arm + "_arc_spec_renderer"],
            width=params[frame_type + "_" + arm + "_arc_spec_width"],
            height=params[frame_type + "_" + arm + "_arc_spec_height"],
            return_jsonstring=params[
                frame_type + "_" + arm + "_arc_spec_return_jsonstring"
            ],
            save_fig=params[frame_type + "_" + arm + "_arc_spec_save_fig"],
            fig_type=params[frame_type + "_" + arm + "_arc_spec_fig_type"],
            filename=os.path.join(
                output_folder_path,
                params[frame_type + "_" + arm + "_arc_spec_filename"]
                + "_"
                + frame_type,
            ),
            open_iframe=params[
                frame_type + "_" + arm + "_arc_spec_open_iframe"
            ],
        )

        # Get the traces for fringe removal
        trace_rectified = copy.deepcopy(twodspec.spectrum_list[0].trace)
        trace_sigma_rectified = copy.deepcopy(
            twodspec.spectrum_list[0].trace_sigma
        )

        if arm == "red":
            # Extract the red flat
            flat = spectral_reduction.TwoDSpec(
                flat_combined_CCDdata.data,
                header=flat_header[0],
                spatial_mask=red_spatial_mask,
                spec_mask=red_spec_mask,
                cosmicray=True,
                readnoise=float(light_temp.header["RDNOISE"]),
                gain=float(light_temp.header["GAIN"]),
                log_file_name=os.path.join(log_file_folder, log_file_name),
            )

            # Add the trace for rectification
            flat.add_trace(trace, trace_sigma)
            flat.get_rectification(coeff=twodspec.rec_coeff)
            flat.apply_rectification()

            # Add the trace for fringe removal
            flat.add_trace(trace_rectified, trace_sigma_rectified)

            # Force extraction from the flat for fringe correction
            flat.ap_extract(apwidth=10, skywidth=0, display=False)

            fringe_correction(twodspec, flat)

    target_name[frame_type] = (
        img[frame_type].light_header[0]["OBJECT"].upper().replace(" ", "")
    )

science_red = extracted_twodspec["science"]["red"]
standard_red = extracted_twodspec["standard"]["red"]
science_blue = extracted_twodspec["science"]["blue"]
standard_blue = extracted_twodspec["standard"]["blue"]

if params["load_standard_target"] is None:
    standard_name = target_name["standard"].lower()
else:
    standard_name = params["load_standard_target"]

science_name = target_name["science"]

#
#
#
# Start handling 1D spectra here
#
#
#
red_onedspec = spectral_reduction.OneDSpec(
    log_file_name=os.path.join(log_file_folder, log_file_name)
)

# Red spectrum first
red_onedspec.from_twodspec(science_red, stype="science")
red_onedspec.from_twodspec(standard_red, stype="standard")

# Find the peaks of the arc
red_onedspec.find_arc_lines(
    prominence=params["red_find_arc_lines_prominence"],
    top_n_peaks=params["red_find_arc_lines_top_n_peaks"],
    distance=params["red_find_arc_lines_distance"],
    refine=params["red_find_arc_lines_refine"],
    refine_window_width=params["red_find_arc_lines_refine_window_width"],
    display=params["red_find_arc_lines_display"],
    width=params["red_find_arc_lines_width"],
    height=params["red_find_arc_lines_height"],
    return_jsonstring=params["red_find_arc_lines_return_jsonstring"],
    renderer=params["red_find_arc_lines_renderer"],
    save_fig=params["red_find_arc_lines_save_fig"],
    fig_type=params["red_find_arc_lines_fig_type"],
    filename=os.path.join(
        output_folder_path, params["red_find_arc_lines_filename"] + "_science"
    ),
    open_iframe=params["red_find_arc_lines_open_iframe"],
    stype="science",
)
red_onedspec.find_arc_lines(
    prominence=params["red_find_arc_lines_prominence"],
    top_n_peaks=params["red_find_arc_lines_top_n_peaks"],
    distance=params["red_find_arc_lines_distance"],
    refine=params["red_find_arc_lines_refine"],
    refine_window_width=params["red_find_arc_lines_refine_window_width"],
    display=params["red_find_arc_lines_display"],
    width=params["red_find_arc_lines_width"],
    height=params["red_find_arc_lines_height"],
    return_jsonstring=params["red_find_arc_lines_return_jsonstring"],
    renderer=params["red_find_arc_lines_renderer"],
    save_fig=params["red_find_arc_lines_save_fig"],
    fig_type=params["red_find_arc_lines_fig_type"],
    filename=os.path.join(
        output_folder_path, params["red_find_arc_lines_filename"] + "_standard"
    ),
    open_iframe=params["red_find_arc_lines_open_iframe"],
    stype="standard",
)

# Configure the wavelength calibrator
red_onedspec.initialise_calibrator(stype="science+standard")

red_onedspec.add_user_atlas(
    elements=element_Hg_red, wavelengths=atlas_Hg_red, stype="science+standard"
)
red_onedspec.add_user_atlas(
    elements=element_Ar_red, wavelengths=atlas_Ar_red, stype="science+standard"
)

red_onedspec.set_hough_properties(
    num_slopes=params["red_set_hough_properties_num_slopes"],
    xbins=params["red_set_hough_properties_xbins"],
    ybins=params["red_set_hough_properties_ybins"],
    min_wavelength=params["red_set_hough_properties_min_wavelength"],
    max_wavelength=params["red_set_hough_properties_max_wavelength"],
    range_tolerance=params["red_set_hough_properties_range_tolerance"],
    linearity_tolerance=params["red_set_hough_properties_linearity_tolerance"],
    stype="science+standard",
)
red_onedspec.set_ransac_properties(
    sample_size=params["red_set_ransac_properties_sample_size"],
    top_n_candidate=params["red_set_ransac_properties_top_n_candidate"],
    linear=params["red_set_ransac_properties_linear"],
    filter_close=params["red_set_ransac_properties_filter_close"],
    ransac_tolerance=params["red_set_ransac_properties_ransac_tolerance"],
    candidate_weighted=params["red_set_ransac_properties_candidate_weighted"],
    hough_weight=params["red_set_ransac_properties_hough_weight"],
    minimum_matches=params["red_set_ransac_properties_minimum_matches"],
    minimum_peak_utilisation=params[
        "red_set_ransac_properties_minimum_peak_utilisation"
    ],
    minimum_fit_error=params["red_set_ransac_properties_minimum_fit_error"],
    stype="science+standard",
)
red_onedspec.do_hough_transform(
    brute_force=params["red_do_hough_transform_brute_force"],
    stype="science+standard",
)

# Solve for the pixel-to-wavelength solution
try:
    red_onedspec.fit(
        max_tries=params["red_fit_max_tries"],
        fit_deg=params["red_fit_fit_deg"],
        fit_coeff=params["red_fit_fit_coeff"],
        fit_tolerance=params["red_fit_fit_tolerance"],
        fit_type=params["red_fit_fit_type"],
        candidate_tolerance=params["red_fit_candidate_tolerance"],
        brute_force=params["red_fit_brute_force"],
        progress=params["red_fit_progress"],
        return_solution=params["red_fit_return_solution"],
        display=params["red_fit_display"],
        renderer=params["red_fit_renderer"],
        save_fig=params["red_fit_save_fig"],
        fig_type=params["red_fit_fig_type"],
        filename=os.path.join(
            output_folder_path, params["red_fit_filename"] + "_science"
        ),
        stype="science",
    )
except:
    try:
        red_onedspec.set_ransac_properties(
            minimum_matches=params["red_set_ransac_properties_minimum_matches"]
            - 1,
            stype="science",
        )
        red_onedspec.fit(
            max_tries=params["red_fit_max_tries"],
            fit_deg=params["red_fit_fit_deg"],
            fit_coeff=params["red_fit_fit_coeff"],
            fit_tolerance=params["red_fit_fit_tolerance"],
            fit_type=params["red_fit_fit_type"],
            candidate_tolerance=params["red_fit_candidate_tolerance"],
            brute_force=params["red_fit_brute_force"],
            progress=params["red_fit_progress"],
            return_solution=params["red_fit_return_solution"],
            display=params["red_fit_display"],
            renderer=params["red_fit_renderer"],
            save_fig=params["red_fit_save_fig"],
            fig_type=params["red_fit_fig_type"],
            filename=os.path.join(
                output_folder_path, params["red_fit_filename"] + "_science"
            ),
            stype="science",
        )
    except:
        if params["hemisphere"] == "north":
            red_onedspec.add_fit_coeff(
                [
                    4768.18,
                    3.42711,
                    3.46385e-6,
                    1.09948e-7,
                    -1.04667e-10,
                    2.96244e-14,
                ],
                stype="science",
            )
        else:
            red_onedspec.add_fit_coeff(
                [
                    4307.37,
                    4.67048,
                    -0.00308638,
                    3.66622e-6,
                    -2.03226e-9,
                    4.27453e-13,
                ],
                stype="science",
            )

try:
    red_onedspec.fit(
        max_tries=params["red_fit_max_tries"],
        fit_deg=params["red_fit_fit_deg"],
        fit_coeff=params["red_fit_fit_coeff"],
        fit_tolerance=params["red_fit_fit_tolerance"],
        fit_type=params["red_fit_fit_type"],
        candidate_tolerance=params["red_fit_candidate_tolerance"],
        brute_force=params["red_fit_brute_force"],
        progress=params["red_fit_progress"],
        return_solution=params["red_fit_return_solution"],
        display=params["red_fit_display"],
        renderer=params["red_fit_renderer"],
        save_fig=params["red_fit_save_fig"],
        fig_type=params["red_fit_fig_type"],
        filename=os.path.join(
            output_folder_path, params["red_fit_filename"] + "_standard"
        ),
        stype="standard",
    )
except:
    try:
        red_onedspec.set_ransac_properties(
            minimum_matches=params["red_set_ransac_properties_minimum_matches"]
            - 1,
            stype="standard",
        )
        red_onedspec.fit(
            max_tries=params["red_fit_max_tries"],
            fit_deg=params["red_fit_fit_deg"],
            fit_coeff=params["red_fit_fit_coeff"],
            fit_tolerance=params["red_fit_fit_tolerance"],
            fit_type=params["red_fit_fit_type"],
            candidate_tolerance=params["red_fit_candidate_tolerance"],
            brute_force=params["red_fit_brute_force"],
            progress=params["red_fit_progress"],
            return_solution=params["red_fit_return_solution"],
            display=params["red_fit_display"],
            renderer=params["red_fit_renderer"],
            save_fig=params["red_fit_save_fig"],
            fig_type=params["red_fit_fig_type"],
            filename=os.path.join(
                output_folder_path, params["red_fit_filename"] + "_standard"
            ),
            stype="standard",
        )
    except:
        if params["hemisphere"] == "north":
            red_onedspec.add_fit_coeff(
                [
                    4768.18,
                    3.42711,
                    3.46385e-6,
                    1.09948e-7,
                    -1.04667e-10,
                    2.96244e-14,
                ],
                stype="standard",
            )
        else:
            red_onedspec.add_fit_coeff(
                [
                    4307.37,
                    4.67048,
                    -0.00308638,
                    3.66622e-6,
                    -2.03226e-9,
                    4.27453e-13,
                ],
                stype="standard",
            )

# Apply the wavelength calibration and display it
red_onedspec.apply_wavelength_calibration(
    stype="science+standard",
)

red_onedspec.load_standard(
    target=standard_name,
    library=params["load_standard_library"],
    ftype=params["load_standard_ftype"],
    cutoff=params["load_standard_cutoff"],
)
red_onedspec.get_sensitivity(
    k=params["red_get_sensitivity_k"],
    method=params["red_get_sensitivity_method"],
    mask_range=params["red_get_sensitivity_mask_range"],
    mask_fit_order=params["red_get_sensitivity_mask_fit_order"],
    mask_fit_size=params["red_get_sensitivity_mask_fit_size"],
    smooth=params["red_get_sensitivity_smooth"],
    slength=params["red_get_sensitivity_slength"],
    sorder=params["red_get_sensitivity_sorder"],
    return_function=params["red_get_sensitivity_return_function"],
    sens_deg=params["red_get_sensitivity_sens_deg"],
)

red_onedspec.inspect_sensitivity(
    display=params["red_inspect_sensitivity_display"],
    renderer=params["red_inspect_sensitivity_renderer"],
    width=params["red_inspect_sensitivity_width"],
    height=params["red_inspect_sensitivity_height"],
    return_jsonstring=params["red_inspect_sensitivity_return_jsonstring"],
    save_fig=params["red_inspect_sensitivity_save_fig"],
    fig_type=params["red_inspect_sensitivity_fig_type"],
    filename=os.path.join(
        output_folder_path, params["red_inspect_sensitivity_filename"]
    ),
)

red_onedspec.apply_flux_calibration(
    inspect=params["red_apply_flux_calibration_inspect"],
    wave_min=params["red_apply_flux_calibration_wave_min"],
    wave_max=params["red_apply_flux_calibration_wave_max"],
    display=params["red_apply_flux_calibration_display"],
    renderer=params["red_apply_flux_calibration_renderer"],
    width=params["red_apply_flux_calibration_width"],
    height=params["red_apply_flux_calibration_height"],
    return_jsonstring=params["red_apply_flux_calibration_return_jsonstring"],
    save_fig=params["red_apply_flux_calibration_save_fig"],
    fig_type=params["red_apply_flux_calibration_fig_type"],
    filename=os.path.join(
        output_folder_path, params["red_apply_flux_calibration_filename"]
    ),
)
red_onedspec.get_telluric_profile(
    mask_range=params["red_get_telluric_profile_mask_range"]
)

red_onedspec.get_telluric_strength(
    factor=params["red_get_telluric_correction_factor"]
)

red_onedspec.inspect_telluric_profile(
    display=params["red_inspect_telluric_profile_display"],
    renderer=params["red_inspect_telluric_profile_renderer"],
    width=params["red_inspect_telluric_profile_width"],
    height=params["red_inspect_telluric_profile_height"],
    return_jsonstring=params["red_inspect_telluric_profile_return_jsonstring"],
    save_fig=params["red_inspect_telluric_profile_save_fig"],
    fig_type=params["red_inspect_telluric_profile_fig_type"],
    filename=os.path.join(
        output_folder_path, params["red_inspect_telluric_profile_filename"]
    ),
)

red_onedspec.inspect_telluric_correction(
    factor=params["red_inspect_telluric_correction_factor"],
    display=params["red_inspect_telluric_correction_display"],
    renderer=params["red_inspect_telluric_correction_renderer"],
    width=params["red_inspect_telluric_correction_width"],
    height=params["red_inspect_telluric_correction_height"],
    return_jsonstring=params[
        "red_inspect_telluric_correction_return_jsonstring"
    ],
    save_fig=params["red_inspect_telluric_correction_save_fig"],
    fig_type=params["red_inspect_telluric_correction_fig_type"],
    filename=os.path.join(
        output_folder_path, params["red_inspect_telluric_correction_filename"]
    ),
)


red_onedspec.apply_telluric_correction(
    factor=params["red_apply_telluric_correction_factor"]
)

# Blue spectrum here
blue_onedspec = spectral_reduction.OneDSpec(
    log_file_name=os.path.join(log_file_folder, log_file_name)
)

# Then the blue spectrum
blue_onedspec.from_twodspec(science_blue, stype="science")
blue_onedspec.from_twodspec(standard_blue, stype="standard")

# Find the peaks of the arc
blue_onedspec.find_arc_lines(
    prominence=params["blue_find_arc_lines_prominence"],
    top_n_peaks=params["blue_find_arc_lines_top_n_peaks"],
    distance=params["blue_find_arc_lines_distance"],
    refine=params["blue_find_arc_lines_refine"],
    refine_window_width=params["blue_find_arc_lines_refine_window_width"],
    display=params["blue_find_arc_lines_display"],
    width=params["blue_find_arc_lines_width"],
    height=params["blue_find_arc_lines_height"],
    return_jsonstring=params["blue_find_arc_lines_return_jsonstring"],
    renderer=params["blue_find_arc_lines_renderer"],
    save_fig=params["blue_find_arc_lines_save_fig"],
    fig_type=params["blue_find_arc_lines_fig_type"],
    filename=os.path.join(
        output_folder_path, params["blue_find_arc_lines_filename"] + "_science"
    ),
    open_iframe=params["blue_find_arc_lines_open_iframe"],
    stype="science",
)
blue_onedspec.find_arc_lines(
    prominence=params["blue_find_arc_lines_prominence"],
    top_n_peaks=params["blue_find_arc_lines_top_n_peaks"],
    distance=params["blue_find_arc_lines_distance"],
    refine=params["blue_find_arc_lines_refine"],
    refine_window_width=params["blue_find_arc_lines_refine_window_width"],
    display=params["blue_find_arc_lines_display"],
    width=params["blue_find_arc_lines_width"],
    height=params["blue_find_arc_lines_height"],
    return_jsonstring=params["blue_find_arc_lines_return_jsonstring"],
    renderer=params["blue_find_arc_lines_renderer"],
    save_fig=params["blue_find_arc_lines_save_fig"],
    fig_type=params["blue_find_arc_lines_fig_type"],
    filename=os.path.join(
        output_folder_path,
        params["blue_find_arc_lines_filename"] + "_standard",
    ),
    open_iframe=params["blue_find_arc_lines_open_iframe"],
    stype="standard",
)


# Configure the wavelength calibrator
blue_onedspec.initialise_calibrator(stype="science+standard")

blue_onedspec.add_user_atlas(
    elements=element_Hg_blue,
    wavelengths=atlas_Hg_blue,
    stype="science+standard",
)
blue_onedspec.add_user_atlas(
    elements=element_Ar_blue,
    wavelengths=atlas_Ar_blue,
    stype="science+standard",
)
blue_onedspec.add_user_atlas(
    elements=element_Zn_blue,
    wavelengths=atlas_Zn_blue,
    stype="science+standard",
)

blue_onedspec.set_hough_properties(
    num_slopes=params["blue_set_hough_properties_num_slopes"],
    xbins=params["blue_set_hough_properties_xbins"],
    ybins=params["blue_set_hough_properties_ybins"],
    min_wavelength=params["blue_set_hough_properties_min_wavelength"],
    max_wavelength=params["blue_set_hough_properties_max_wavelength"],
    range_tolerance=params["blue_set_hough_properties_range_tolerance"],
    linearity_tolerance=params[
        "blue_set_hough_properties_linearity_tolerance"
    ],
    stype="science+standard",
)
blue_onedspec.set_ransac_properties(
    sample_size=params["blue_set_ransac_properties_sample_size"],
    top_n_candidate=params["blue_set_ransac_properties_top_n_candidate"],
    linear=params["blue_set_ransac_properties_linear"],
    filter_close=params["blue_set_ransac_properties_filter_close"],
    ransac_tolerance=params["blue_set_ransac_properties_ransac_tolerance"],
    candidate_weighted=params["blue_set_ransac_properties_candidate_weighted"],
    hough_weight=params["blue_set_ransac_properties_hough_weight"],
    minimum_matches=params["blue_set_ransac_properties_minimum_matches"],
    minimum_peak_utilisation=params[
        "blue_set_ransac_properties_minimum_peak_utilisation"
    ],
    minimum_fit_error=params["blue_set_ransac_properties_minimum_fit_error"],
    stype="science+standard",
)
blue_onedspec.do_hough_transform(
    brute_force=params["blue_do_hough_transform_brute_force"],
    stype="science+standard",
)

# Solve for the pixel-to-wavelength solution
try:
    blue_onedspec.fit(
        max_tries=params["blue_fit_max_tries"],
        fit_deg=params["blue_fit_fit_deg"],
        fit_coeff=params["blue_fit_fit_coeff"],
        fit_tolerance=params["blue_fit_fit_tolerance"],
        fit_type=params["blue_fit_fit_type"],
        candidate_tolerance=params["blue_fit_candidate_tolerance"],
        brute_force=params["blue_fit_brute_force"],
        progress=params["blue_fit_progress"],
        return_solution=params["blue_fit_return_solution"],
        display=params["blue_fit_display"],
        renderer=params["blue_fit_renderer"],
        save_fig=params["blue_fit_save_fig"],
        fig_type=params["blue_fit_fig_type"],
        filename=os.path.join(
            output_folder_path, params["blue_fit_filename"] + "_science"
        ),
        stype="science",
    )
except:
    if params["hemisphere"] == "north":
        blue_onedspec.add_fit_coeff(
            [3323.6, 1.71555, -5.14076e-5, 7.39966e18, -2.24418e-11],
            stype="science",
        )
    else:
        blue_onedspec.add_fit_coeff(
            [3176.54, 1.76596, -0.000130639, 1.21175e-7, -3.33912e-11],
            stype="science",
        )

try:
    blue_onedspec.fit(
        max_tries=params["blue_fit_max_tries"],
        fit_deg=params["blue_fit_fit_deg"],
        fit_coeff=params["blue_fit_fit_coeff"],
        fit_tolerance=params["blue_fit_fit_tolerance"],
        fit_type=params["blue_fit_fit_type"],
        candidate_tolerance=params["blue_fit_candidate_tolerance"],
        brute_force=params["blue_fit_brute_force"],
        progress=params["blue_fit_progress"],
        return_solution=params["blue_fit_return_solution"],
        display=params["blue_fit_display"],
        renderer=params["blue_fit_renderer"],
        save_fig=params["blue_fit_save_fig"],
        fig_type=params["blue_fit_fig_type"],
        filename=os.path.join(
            output_folder_path, params["blue_fit_filename"] + "_standard"
        ),
        stype="standard",
    )
except:
    if params["hemisphere"] == "north":
        blue_onedspec.add_fit_coeff(
            [3323.6, 1.71555, -5.14076e-5, 7.39966e18, -2.24418e-11],
            stype="standard",
        )
    else:
        blue_onedspec.add_fit_coeff(
            [3176.54, 1.76596, -0.000130639, 1.21175e-7, -3.33912e-11],
            stype="standard",
        )

# Apply the wavelength calibration and display it
blue_onedspec.apply_wavelength_calibration(
    stype="science+standard",
)

blue_onedspec.load_standard(
    target=standard_name,
    library=params["load_standard_library"],
    ftype=params["load_standard_ftype"],
    cutoff=params["load_standard_cutoff"],
)
blue_onedspec.get_sensitivity(
    k=params["blue_get_sensitivity_k"],
    method=params["blue_get_sensitivity_method"],
    mask_range=params["blue_get_sensitivity_mask_range"],
    mask_fit_order=params["blue_get_sensitivity_mask_fit_order"],
    mask_fit_size=params["blue_get_sensitivity_mask_fit_size"],
    smooth=params["blue_get_sensitivity_smooth"],
    slength=params["blue_get_sensitivity_slength"],
    sorder=params["blue_get_sensitivity_sorder"],
    return_function=params["blue_get_sensitivity_return_function"],
    sens_deg=params["blue_get_sensitivity_sens_deg"],
)

blue_onedspec.inspect_sensitivity(
    display=params["blue_inspect_sensitivity_display"],
    renderer=params["blue_inspect_sensitivity_renderer"],
    width=params["blue_inspect_sensitivity_width"],
    height=params["blue_inspect_sensitivity_height"],
    return_jsonstring=params["blue_inspect_sensitivity_return_jsonstring"],
    save_fig=params["blue_inspect_sensitivity_save_fig"],
    fig_type=params["blue_inspect_sensitivity_fig_type"],
    filename=os.path.join(
        output_folder_path, params["blue_inspect_sensitivity_filename"]
    ),
)

blue_onedspec.apply_flux_calibration(
    inspect=params["blue_apply_flux_calibration_inspect"],
    wave_min=params["blue_apply_flux_calibration_wave_min"],
    wave_max=params["blue_apply_flux_calibration_wave_max"],
    display=params["blue_apply_flux_calibration_display"],
    renderer=params["blue_apply_flux_calibration_renderer"],
    width=params["blue_apply_flux_calibration_width"],
    height=params["blue_apply_flux_calibration_height"],
    return_jsonstring=params["blue_apply_flux_calibration_return_jsonstring"],
    save_fig=params["blue_apply_flux_calibration_save_fig"],
    fig_type=params["blue_apply_flux_calibration_fig_type"],
    filename=os.path.join(
        output_folder_path, params["blue_apply_flux_calibration_filename"]
    ),
)


#
#
#
# second pass to correct for the small scale detector response
#
#
#
red_std_spec = red_onedspec.standard_spectrum_list[0]
red_sci_spec = red_onedspec.science_spectrum_list[0]

red_literature_flux = spectres(
    red_std_spec.wave,
    red_std_spec.wave_literature,
    red_std_spec.flux_literature,
)

red_second_pass_correction = red_std_spec.flux / red_literature_flux
red_second_pass_correction = medfilt(red_second_pass_correction, 5)

red_std_spec.flux /= red_second_pass_correction
red_std_spec.flux_err /= red_second_pass_correction
red_std_spec.flux_sky /= red_second_pass_correction

red_sci_spec.flux /= red_second_pass_correction
red_sci_spec.flux_err /= red_second_pass_correction
red_sci_spec.flux_sky /= red_second_pass_correction

blue_std_spec = blue_onedspec.standard_spectrum_list[0]
blue_sci_spec = blue_onedspec.science_spectrum_list[0]

blue_literature_flux = spectres(
    blue_std_spec.wave,
    blue_std_spec.wave_literature,
    blue_std_spec.flux_literature,
)

blue_second_pass_correction = blue_std_spec.flux / blue_literature_flux
blue_second_pass_correction = medfilt(blue_second_pass_correction, 15)

blue_std_spec.flux /= blue_second_pass_correction
blue_std_spec.flux_err /= blue_second_pass_correction
blue_std_spec.flux_sky /= blue_second_pass_correction

blue_sci_spec.flux /= blue_second_pass_correction
blue_sci_spec.flux_err /= blue_second_pass_correction
blue_sci_spec.flux_sky /= blue_second_pass_correction

#
#
#
# Atmospheric extinction correction
#
#
#
red_onedspec.set_atmospheric_extinction(
    location=params["red_set_atmospheric_extinction_location"],
    kind=params["red_set_atmospheric_extinction_kind"],
    fill_value=params["red_set_atmospheric_extinction_fill_value"],
)

blue_onedspec.set_atmospheric_extinction(
    location=params["blue_set_atmospheric_extinction_location"],
    kind=params["blue_set_atmospheric_extinction_kind"],
    fill_value=params["blue_set_atmospheric_extinction_fill_value"],
)

red_onedspec.apply_atmospheric_extinction_correction(
    science_airmass=img["science"].light_header[0]["AIRMASS"],
    standard_airmass=img["standard"].light_header[0]["AIRMASS"],
)

blue_onedspec.apply_atmospheric_extinction_correction(
    science_airmass=img["science"].light_header[0]["AIRMASS"],
    standard_airmass=img["standard"].light_header[0]["AIRMASS"],
)

#
#
#
# Save figures
#
#
#
blue_science_filename = params[
    "blue_inspect_reduced_spectrum_science_filename"
]
if blue_science_filename is None:
    blue_science_filename = os.path.join(
        output_folder_path, science_name + "_blue"
    )
blue_onedspec.inspect_reduced_spectrum(
    wave_min=params["blue_inspect_reduced_spectrum_wave_min"],
    wave_max=params["blue_inspect_reduced_spectrum_wave_max"],
    display=params["blue_inspect_reduced_spectrum_display"],
    renderer=params["blue_inspect_reduced_spectrum_renderer"],
    width=params["blue_inspect_reduced_spectrum_width"],
    height=params["blue_inspect_reduced_spectrum_height"],
    save_fig=params["blue_inspect_reduced_spectrum_save_fig"],
    fig_type=params["blue_inspect_reduced_spectrum_fig_type"],
    filename=os.path.join(output_folder_path, blue_science_filename),
    stype="science",
)

red_science_filename = params["red_inspect_reduced_spectrum_science_filename"]
if red_science_filename is None:
    red_science_filename = os.path.join(
        output_folder_path, science_name + "_blue"
    )
red_onedspec.inspect_reduced_spectrum(
    wave_min=params["red_inspect_reduced_spectrum_wave_min"],
    wave_max=params["red_inspect_reduced_spectrum_wave_max"],
    display=params["red_inspect_reduced_spectrum_display"],
    renderer=params["red_inspect_reduced_spectrum_renderer"],
    width=params["red_inspect_reduced_spectrum_width"],
    height=params["red_inspect_reduced_spectrum_height"],
    save_fig=params["red_inspect_reduced_spectrum_save_fig"],
    fig_type=params["red_inspect_reduced_spectrum_fig_type"],
    filename=os.path.join(output_folder_path, red_science_filename),
    stype="science",
)

blue_standard_filename = params[
    "blue_inspect_reduced_spectrum_standard_filename"
]
if blue_standard_filename is None:
    blue_standard_filename = os.path.join(
        output_folder_path, standard_name + "_blue"
    )
blue_onedspec.inspect_reduced_spectrum(
    wave_min=params["blue_inspect_reduced_spectrum_wave_min"],
    wave_max=params["blue_inspect_reduced_spectrum_wave_max"],
    display=params["blue_inspect_reduced_spectrum_display"],
    renderer=params["blue_inspect_reduced_spectrum_renderer"],
    width=params["blue_inspect_reduced_spectrum_width"],
    height=params["blue_inspect_reduced_spectrum_height"],
    save_fig=params["blue_inspect_reduced_spectrum_save_fig"],
    fig_type=params["blue_inspect_reduced_spectrum_fig_type"],
    filename=os.path.join(output_folder_path, blue_standard_filename),
    stype="standard",
)

red_standard_filename = params[
    "red_inspect_reduced_spectrum_standard_filename"
]
if red_standard_filename is None:
    red_standard_filename = os.path.join(
        output_folder_path, standard_name + "_blue"
    )
red_onedspec.inspect_reduced_spectrum(
    wave_min=params["red_inspect_reduced_spectrum_wave_min"],
    wave_max=params["red_inspect_reduced_spectrum_wave_max"],
    display=params["red_inspect_reduced_spectrum_display"],
    renderer=params["red_inspect_reduced_spectrum_renderer"],
    width=params["red_inspect_reduced_spectrum_width"],
    height=params["red_inspect_reduced_spectrum_height"],
    save_fig=params["red_inspect_reduced_spectrum_save_fig"],
    fig_type=params["red_inspect_reduced_spectrum_fig_type"],
    filename=os.path.join(output_folder_path, red_standard_filename),
    stype="standard",
)


#
#
#
# Resampling before saving
#
#
#
red_onedspec.resample(
    wave_start=params["red_apply_wavelength_calibration_wave_start"],
    wave_end=params["red_apply_wavelength_calibration_wave_end"],
    wave_bin=params["red_apply_wavelength_calibration_wave_bin"],
)
blue_onedspec.resample(
    wave_start=params["blue_apply_wavelength_calibration_wave_start"],
    wave_end=params["blue_apply_wavelength_calibration_wave_end"],
    wave_bin=params["blue_apply_wavelength_calibration_wave_bin"],
)


#
#
#
# Saving files here
#
#
#
if params["red_create_fits"]:
    red_onedspec.create_fits(
        output=params["red_create_fits_output"],
        recreate=params["red_create_fits_recreate"],
        empty_primary_hdu=params["red_create_fits_empty_primary_hdu"],
    )

if params["blue_create_fits"]:
    blue_onedspec.create_fits(
        output=params["blue_create_fits_output"],
        recreate=params["blue_create_fits_recreate"],
        empty_primary_hdu=params["blue_create_fits_empty_primary_hdu"],
    )

if params["red_save_fits"]:
    red_onedspec.save_fits(
        output=params["red_save_fits_output"],
        filename=os.path.join(
            output_folder_path, params["red_save_fits_filename"]
        ),
        recreate=params["red_save_fits_recreate"],
        empty_primary_hdu=params["red_save_fits_empty_primary_hdu"],
        overwrite=params["red_save_fits_overwrite"],
    )

if params["blue_save_fits"]:
    blue_onedspec.save_fits(
        output=params["blue_save_fits_output"],
        filename=os.path.join(
            output_folder_path, params["blue_save_fits_filename"]
        ),
        recreate=params["blue_save_fits_recreate"],
        empty_primary_hdu=params["blue_save_fits_empty_primary_hdu"],
        overwrite=params["blue_save_fits_overwrite"],
    )

if params["red_save_csv"]:
    red_onedspec.save_csv(
        output=params["red_save_csv_output"],
        filename=os.path.join(
            output_folder_path, params["red_save_csv_filename"]
        ),
        recreate=params["red_save_csv_recreate"],
        overwrite=params["red_save_csv_overwrite"],
    )

if params["blue_save_csv"]:
    blue_onedspec.save_csv(
        output=params["blue_save_csv_output"],
        filename=os.path.join(
            output_folder_path, params["blue_save_csv_filename"]
        ),
        recreate=params["blue_save_csv_recreate"],
        overwrite=params["blue_save_csv_overwrite"],
    )

#
#
#
# Stitching up the red and blue and save csv and png
#
#
#
wave_red = red_onedspec.science_spectrum_list[0].wave_resampled

flux_red = red_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_telluric_corrected
flux_red_err = red_onedspec.science_spectrum_list[
    0
].flux_err_resampled_atm_ext_telluric_corrected

# blue is NOT telluric corrected as it is not needed
wave_blue = blue_onedspec.science_spectrum_list[0].wave_resampled
flux_blue = blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected
flux_blue_err = blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected

# trim the last few hundred A from the blue and the first few hundred A from
# the red in the combined spectrum
red_limit = 5000
# This value has to match the
blue_limit = 5750

blue_mask = (wave_blue >= red_limit) & (wave_blue <= blue_limit)
red_mask = (wave_red >= red_limit) & (wave_red <= blue_limit)

# resample the red to match blue resolution
flux_red_resampled, flux_red_resampled_err = spectres(
    wave_blue[blue_mask],
    wave_red[red_mask],
    flux_red[red_mask],
    flux_red_err[red_mask],
)

flux_weighted_combine = (
    flux_red_resampled / flux_red_resampled_err**2.0
    + flux_blue[blue_mask] / flux_blue_err[blue_mask] ** 2.0
) / (
    1.0 / flux_red_resampled_err**2.0 + 1.0 / flux_blue_err[blue_mask] ** 2.0
)

flux_error_weighted_combine = 1.0 / np.sqrt(
    1.0 / flux_red_resampled_err**2.0 + 1.0 / flux_blue_err[blue_mask] ** 2.0
)

flux_combined = np.concatenate(
    (
        flux_blue[wave_blue < red_limit],
        flux_weighted_combine,
        flux_red[wave_red > blue_limit],
    )
)
wave_combined = np.concatenate(
    (wave_blue[wave_blue < blue_limit], wave_red[wave_red > blue_limit])
)

flux_error_combined = np.concatenate(
    (
        flux_blue_err[wave_blue < red_limit],
        flux_error_weighted_combine,
        flux_red_err[wave_red > blue_limit],
    )
)

# https://lco.global/observatory/instruments/floyds/
# pixel scale is 1.74
wave_out = np.arange(3300.0, 10000.0, 1.74)
flux_out, flux_error_out = spectres(
    wave_out, wave_combined, flux_combined, flux_error_combined
)

plt.figure(1, figsize=(16, 8))
plt.clf()
plt.plot(wave_blue, flux_blue, color="blue", label="Blue arm (ASPIRED)")
plt.plot(wave_red, flux_red, color="red", label="Red arm (ASPIRED)")
plt.plot(wave_out, flux_out, color="black", label="Weighted Combined")
plt.fill_between(
    wave_out,
    flux_out - flux_error_out,
    flux_out + flux_error_out,
    color="gray",
    label=r"Weighted Combined Noise (1$\sigma$)",
    zorder=1,
)
plt.xlim(3300.0, 10000.0)
plt.ylim(
    np.nanpercentile(
        flux_combined[wave_combined < 9000][
            np.isfinite(flux_combined[wave_combined < 9000])
        ],
        0.25,
    ),
    np.nanpercentile(
        flux_combined[wave_combined < 9000][
            np.isfinite(flux_combined[wave_combined < 9000])
        ],
        99.5,
    ),
)
plt.xlabel("Wavelength / A")
plt.ylabel("Flux / (erg / s / cm / cm / A)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.title(science_name)
plt.savefig(
    os.path.join(
        output_folder_path,
        science_name.replace(" ", "_")
        + "_both_arms_{}.png".format(params["output_file_name_suffix"]),
    )
)
np.savetxt(
    os.path.join(
        output_folder_path,
        science_name.replace(" ", "_")
        + "_{}.csv".format(params["output_file_name_suffix"]),
    ),
    np.column_stack((wave_combined, flux_combined, flux_error_combined)),
    delimiter=",",
)

utc_time = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
#
#
#
# Save fits with header suitable for SNEx
# https://lco.global/observatory/instruments/floyds/
# resolution is ~1.74 A / pix
#
blue_onedspec.create_fits(output="flux_resampled_atm_ext_corrected")

blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].data = flux_out

blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header = img[
    "science"
].light_header[
    0
]

blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["nstack"] = len(
    light_path
)

for i in range(len(light_path)):
    blue_onedspec.science_spectrum_list[
        0
    ].flux_resampled_atm_ext_corrected_hdulist[0].header[
        "frame{}".format(i + 1)
    ] = light_path[
        i
    ].split(
        os.path.sep
    )[
        -1
    ]

blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header[
    "SLIT"
] = light_temp.header[
    "APERWID"
]
blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header[
    "EXPTIME"
] = total_exposure_time

blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["CTYPE1"] = "Wavelength"
blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["CUNIT1"] = "Angstroms"
blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["CRVAL1"] = 3.300e03
blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["CDELT1"] = 1.74e00
blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["CD1_1"] = 1.74e00
blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["CRPIX1"] = 1.00e00
blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["CTYPE2"] = "LINEAR"
blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["CUNIT2"] = "Pixels  "
blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["CD2_2"] = 0.00e00

blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["REDUCER"] = params[
    "reducer"
]
blue_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_corrected_hdulist[0].header["REDUCET"] = utc_time

output = fits.PrimaryHDU(
    blue_onedspec.science_spectrum_list[0]
    .flux_resampled_atm_ext_corrected_hdulist[0]
    .data,
    blue_onedspec.science_spectrum_list[0]
    .flux_resampled_atm_ext_corrected_hdulist[0]
    .header,
)
output.writeto(
    os.path.join(
        output_folder_path,
        science_name.replace(" ", "_")
        + "_{}.fits".format(params["output_file_name_suffix"]),
    ),
    overwrite=True,
)
