#!/usr/bin/python3
"""
This is a helper sciprt to allow (1) getting TNS info, (2) downloading
FLOYDS L1 data from LCO archive, and (3) directly using local files
for reduction

"""

from argparse import ArgumentParser
from collections import OrderedDict
from datetime import timedelta, datetime
import getpass
import glob
import json
import os
import sys

from astropy.io import fits
import numpy as np
import requests
import ruamel.yaml

from query_lco_archive import get_metadata, download_frame
from query_tns import tns_search, get_tns

# Configure the parser
parser = ArgumentParser(
    description="Downloads and Reduce data from archive.lco.global."
)

# Note that the latter half of these arguments are ignored when --local is
# True, which is the default
# When processing locally, ONLY ONE STANDARD can be present in the folder,
# all the science frames will get processed with that standard.
parser.add_argument(
    "--reducer",
    default="Automated",
    help="The reducer's name to be added to the header.",
)
parser.add_argument(
    "--observer",
    default="Automated",
    help="The observer's name to be added to the header.",
)
parser.add_argument(
    "--target_name",
    default=None,
    help="The target name to be queried on the TNS.",
)
parser.add_argument(
    "--ra",
    default=None,
    help="Right Ascension in decimal. Only used if target_name is None.",
)
parser.add_argument(
    "--dec",
    default=None,
    help="Declination in decimal. Only used if target_name is None.",
)
parser.add_argument(
    "--input_folder",
    default=None,
    help="Path to store the raw data products.",
)
parser.add_argument(
    "--output_folder",
    default=None,
    help="Path to store the reduced data products.",
)
parser.add_argument(
    "--local",
    action="store_true",
    help="Set to reduce local data without connecting to TNS/LCO",
)
# If local is True, all the following arguments will get ignored
parser.add_argument("--login", default=None, help="Path to the login details.")
parser.add_argument(
    "--lco_token",
    default=None,
    help=(
        "LCO token. Only used if --login is None. Will ask for one if"
        "needed and nothing is provided."
    ),
)
parser.add_argument(
    "--tns_bot_id",
    default=None,
    help=(
        "TNS Bot ID. Only used if --login is None. Will ask for one if"
        "needed and nothing is provided."
    ),
)
parser.add_argument(
    "--tns_bot_name",
    default=None,
    help=(
        "TNS Bot name. Only used if --login is None. Will ask for one if"
        "needed and nothing is provided."
    ),
)
parser.add_argument(
    "--tns_token",
    default=None,
    help=(
        "TNS token. Only used if --login is None. Will ask for one if"
        "needed and nothing is provided."
    ),
)
parser.add_argument(
    "--date_start",
    default="1900-01-01",
    help="The date of the beginning of the night of the observation.",
)
parser.add_argument(
    "--date_end",
    default="2100-12-31",
    help="The date of the beginning of the night of the observation.",
)
parser.add_argument(
    "--most_recent_only",
    action="store_true",
    help="Set to reduce the most recently collected spectrum only.",
)
args = parser.parse_args()


HERE = os.path.dirname(os.path.realpath(__file__))

# Modify the yaml initialisation
yaml = ruamel.yaml.YAML()
yaml.preserve_quotes = True
yaml_template = os.path.join(HERE, "floyds_template.yaml")


# Load the login information: TNS & SNEx
if args.login is not None:
    with open(args.login, encoding="ascii") as f:
        login_yaml = yaml.load(f)

# choose directory to store the output
input_folder = args.input_folder
output_folder = args.output_folder
target_name = args.target_name
reducer = args.reducer
observer = args.observer

# Handle the 4 combinations of how input_folder and target_name are provided
# Case 1: neither
# - raise error
if (input_folder is None) & (target_name is None):
    raise ValueError("Either input_folder or target_name has to be provided")
# Case 2: only target_name is provided
# - set input_folder to target_name
# - get (ra, dec)
elif input_folder is None:
    input_folder = target_name
    """
    if not args.local:
        if (args.ra is None) or (args.dec is None):
            raise ValueError(
                "Either target_name or (ra & dec) has to be provided."
            )
    """
    if args.ra is None:
        ra = None
    else:
        ra = float(args.ra)
    if args.ra is None:
        dec = None
    else:
        dec = float(args.dec)
# Case 3: only input_folder is provided
# - set target_name to input_name
# - get (ra, dec)
elif target_name is None:
    target_name = input_folder
    """
    if not args.local:
        if (args.ra is None) or (args.dec is None):
            raise ValueError(
                "Either target_name or (ra & dec) has to be provided."
            )
    """
    if args.ra is None:
        ra = None
    else:
        ra = float(args.ra)
    if args.ra is None:
        dec = None
    else:
        dec = float(args.dec)
# Case 4 - both are provided
else:
    if args.local:
        if args.ra is None:
            ra = None
        else:
            ra = float(args.ra)
        if args.ra is None:
            dec = None
        else:
            dec = float(args.dec)
    else:
        print("target name is provided, input ra and dec are ignored.")

        # get the target position if only a name is provided
        TNS = "www.wis-tns.org"
        url_tns_api = "https://" + TNS + "/api/get"

        # TNS token: needed for name resolving to (ra, dec)
        if args.login is not None:
            TNS_BOT_ID = login_yaml["TNS_BOT_ID"]
            TNS_BOT_NAME = login_yaml["TNS_BOT_NAME"]
            TNS_API_KEY = login_yaml["TNS_API_KEY"]
        else:
            TNS_BOT_ID = args.tns_bot_id
            TNS_BOT_NAME = args.tns_bot_name
            TNS_API_KEY = args.tns_token

        # Prompt to ask for the IDs if some are still None.
        if TNS_BOT_ID is None:
            TNS_BOT_ID = getpass.getpass(prompt="Enter TNS_BOT_ID: ")
        if TNS_BOT_NAME is None:
            TNS_BOT_NAME = getpass.getpass(prompt="Enter TNS_BOT_NAME: ")
        if TNS_API_KEY is None:
            TNS_API_KEY = getpass.getpass(prompt="Enter TNS_API_KEY: ")

        search_obj = [("objname", target_name)]
        tns_search_response = tns_search(search_obj)

        # assuming only 1 result...
        get_obj = [
            ("objid", tns_search_response.json()["data"]["reply"][0]["objid"])
        ]
        get_obj_response = get_tns(get_obj)

        ra = get_obj_response.json()["data"]["reply"]["radeg"]
        dec = get_obj_response.json()["data"]["reply"]["decdeg"]


# Get the absolute path to the input folder
if os.path.isabs(input_folder):
    input_folder_abs_path = input_folder
else:
    input_folder_abs_path = os.path.join(HERE, input_folder)

# Get the absolute path to the output folder
# default to be in a folder in the same directory of the input_folder
# with the folder name of the input_folder appeneded with output
if output_folder is None:
    output_folder = input_folder + "_output"

if os.path.isabs(output_folder):
    output_folder_abs_path = output_folder
else:
    output_folder_abs_path = os.path.join(HERE, output_folder)

# if --local is False
# Download the data
if not args.local:
    # get the token
    # LCO token: needed for private data
    if args.login is not None:
        lco_token = login_yaml["lco_token"]

    else:
        lco_token = args.lco_token

    if lco_token == "public":
        authtoken = {}

    else:
        if lco_token is None:
            # Prompt to ask for the token is still None.
            lco_token = getpass.getpass(prompt="Enter LCO API token: ")
            authtoken = {"Authorization": "Token {}".format(lco_token)}

        else:
            authtoken = {"Authorization": "Token {}".format(lco_token)}

    # get the filelist
    # en06 is FLOYDS North at the Haleakala Observatory (OGG)
    # en12 is FLOYDS South at the Siding Spring Observatory (COJ)
    # the search has to be done seprately
    date_start = args.date_start
    date_end = args.date_end

    science_metadata = []
    standard_metadata = []
    request_id = []
    obstype = []
    # The "day" of the beginning of the observation
    day_obs = []
    # The exact time of the start of the observation
    date_obs = []
    instrume_list = []

    # first north (en06) then south (en12)
    for instrume in ["en06", "en12"]:
        # Get the metadata
        _science_metadata = get_metadata(
            authtoken=authtoken,
            limit=1000,
            INSTRUME=instrume,
            start=date_start,
            end=date_end,
            RLEVEL=0,
            covers="POINT({} {})".format(ra, dec),
        )

        # There is only one spectrum from a single group of request_id
        for metadata in _science_metadata:
            request_id.append(metadata["request_id"])
            obstype.append(metadata["OBSTYPE"])
            day_obs.append(metadata["DAY_OBS"])
            date_obs.append(metadata["DATE_OBS"][:-1])
            instrume_list.append(instrume)
        science_metadata += _science_metadata

    request_id = np.array(request_id)
    obstype = np.array(obstype)
    day_obs = np.array(day_obs)
    date_obs = np.array(date_obs)
    instrume_list = np.array(instrume_list)

    # If only requesting the newest data, remove the unwanted results
    if args.most_recent_only:
        on_target = obstype == "SPECTRUM"
        _keeping = np.argmax(request_id[on_target])
        request_id = [np.asarray(request_id[on_target])[_keeping]]
        obstype = [np.asarray(obstype[on_target])[_keeping]]
        day_obs = [np.asarray(day_obs[on_target])[_keeping]]
        date_obs = [np.asarray(date_obs[on_target])[_keeping]]

    # Get the request_id that contains useful spectral data
    request_id_science = [
        i for i, j in zip(request_id, obstype) if j == "SPECTRUM"
    ]

    science_metadata = [
        v for v in science_metadata if v["request_id"] in request_id_science
    ]

    # Look for the standard frames taken closest to the time of observation
    for date, instrume, day in [
        (i, j, l)
        for i, j, k, l in zip(date_obs, instrume_list, day_obs, obstype)
        if l == "SPECTRUM"
    ]:
        _date = datetime.fromisoformat(date)
        day_range = 0.5
        _standard_metadata = []
        while _standard_metadata == []:
            # Get the standard star
            _standard_metadata = get_metadata(
                authtoken=authtoken,
                limit=1000,
                INSTRUME=instrume,
                start=(_date - timedelta(days=day_range)).isoformat(),
                end=(_date + timedelta(days=day_range)).isoformat(),
                PROPID="FLOYDS standards",
                OBSTYPE="SPECTRUM",
                RLEVEL=0,
            )
            day_range += 1
        for i in _standard_metadata:
            print(i["DAY_OBS"])
        # If there are more than 1 frame returned, check for matching day_obs
        # first and then check for the one with the least time difference
        _standard_metadata_day_obs = []
        _time = []
        if len(_standard_metadata) > 1:
            # first check for the day_obs
            for meta in _standard_metadata:
                print(meta["DAY_OBS"])
                if meta["DAY_OBS"] == day_obs:
                    _standard_metadata_day_obs.append(meta)
                if len(_standard_metadata_day_obs) > 1:
                    # Get all the timestamps
                    for meta in _standard_metadata_day_obs:
                        _time.append(
                            datetime.fromisoformat(
                                meta["observation_date"][:-1]
                            )
                        )
                    _time_diff = [i - _date for i in _time]
                    min_time_idx = np.argmin(np.abs(np.array(_time_diff)))
                elif len(_standard_metadata_day_obs) == 1:
                    min_time_idx = 0
                    _time.append(
                        datetime.fromisoformat(meta["observation_date"][:-1])
                    )
                else:
                    # Get all the timestamps
                    for meta in _standard_metadata:
                        print(meta["DAY_OBS"])
                        _time.append(
                            datetime.fromisoformat(
                                meta["observation_date"][:-1]
                            )
                        )
                    _time_diff = [i - _date for i in _time]
                    min_time_idx = np.argmin(np.abs(np.array(_time_diff)))
        else:
            min_time_idx = 0
            _time.append(
                datetime.fromisoformat(
                    _standard_metadata[0]["observation_date"][:-1]
                )
            )

    standard_spectrum_metadata = get_metadata(
        authtoken=authtoken,
        limit=1000,
        INSTRUME=instrume,
        start=(_time[min_time_idx] - timedelta(minutes=30)).isoformat(),
        end=(_time[min_time_idx] + timedelta(minutes=30)).isoformat(),
        PROPID="FLOYDS standards",
        OBSTYPE="SPECTRUM",
        RLEVEL=0,
    )
    if len(standard_spectrum_metadata) > 1:
        standard_spectrum_metadata = [d for d in standard_spectrum_metadata][0]

    # make sure there is an arc...
    standard_arc_metadata = []
    day_range = 0.0
    while standard_arc_metadata == []:
        standard_arc_metadata = get_metadata(
            authtoken=authtoken,
            limit=1000,
            INSTRUME=instrume,
            start=(
                _time[min_time_idx]
                - timedelta(minutes=30)
                - timedelta(days=day_range)
            ).isoformat(),
            end=(
                _time[min_time_idx]
                + timedelta(minutes=30)
                + timedelta(days=day_range)
            ).isoformat(),
            PROPID="FLOYDS standards",
            OBSTYPE="ARC",
            RLEVEL=0,
        )
        day_range += 1
    # make sure there is a flat...
    standard_flat_metadata = []
    day_range = 0.0
    while standard_flat_metadata == []:
        standard_flat_metadata = get_metadata(
            authtoken=authtoken,
            limit=1000,
            INSTRUME=instrume,
            start=(
                _time[min_time_idx]
                - timedelta(minutes=30)
                - timedelta(days=day_range)
            ).isoformat(),
            end=(
                _time[min_time_idx]
                + timedelta(minutes=30)
                + timedelta(days=day_range)
            ).isoformat(),
            PROPID="FLOYDS standards",
            OBSTYPE="SPECTRUM",
            RLEVEL=0,
        )
        if len(standard_spectrum_metadata) > 1:
            standard_spectrum_metadata = [
                d for d in standard_spectrum_metadata
            ][0]

        # make sure there is an arc...
        standard_arc_metadata = []
        day_range = 0.0
        while standard_arc_metadata == []:
            standard_arc_metadata = get_metadata(
                authtoken=authtoken,
                limit=1000,
                INSTRUME=instrume,
                start=(
                    _time[min_time_idx]
                    - timedelta(minutes=30)
                    - timedelta(days=day_range)
                ).isoformat(),
                end=(
                    _time[min_time_idx]
                    + timedelta(minutes=30)
                    + timedelta(days=day_range)
                ).isoformat(),
                PROPID="FLOYDS standards",
                OBSTYPE="ARC",
                RLEVEL=0,
            )
            day_range += 1
        # make sure there is a flat...
        standard_flat_metadata = []
        day_range = 0.0
        while standard_flat_metadata == []:
            standard_flat_metadata = get_metadata(
                authtoken=authtoken,
                limit=1000,
                INSTRUME=instrume,
                start=(
                    _time[min_time_idx] - timedelta(minutes=30)
                ).isoformat(),
                end=(_time[min_time_idx] + timedelta(minutes=30)).isoformat(),
                PROPID="FLOYDS standards",
                OBSTYPE="LAMPFLAT",
                RLEVEL=0,
            )
            day_range += 1
        if len(standard_spectrum_metadata) > 1:
            standard_spectrum_metadata = standard_spectrum_metadata[0]
        if len(standard_flat_metadata) > 1:
            standard_flat_metadata = standard_flat_metadata[0]
        if len(standard_arc_metadata) > 1:
            standard_arc_metadata = standard_arc_metadata[0]
        standard_metadata += standard_spectrum_metadata
        standard_metadata += standard_arc_metadata
        standard_metadata += standard_flat_metadata

    # Pack the light, flat & arc fro science and standard as a dictionary item
    target_list = {}

    # Download the data and distinguish north-south for the rectification
    science_filepath_list = []
    science_filename_list = []
    science_hemisphere_list = []
    standard_filepath_list = []
    standard_filename_list = []
    standard_hemisphere_list = []

    # Download the science and the best matched standard frames
    for science_frame, standard_frame in zip(
        science_metadata, standard_metadata
    ):
        if (
            (".tar.gz" not in science_frame["filename"])
            or ("fits.fz" not in science_frame["filename"])
        ) & (science_frame["request_id"] in request_id_science):
            [
                x.append(y)
                for x, y in zip(
                    [science_filepath_list, science_filename_list],
                    download_frame(
                        frame=science_frame,
                        base_directory=input_folder,
                        no_date=False,
                    ),
                )
            ]
            # Arrange into the target list
            if science_frame["request_id"] not in target_list.keys():
                target_list[science_frame["request_id"]] = {
                    "science": {},
                    "standard": {},
                }
            if science_frame["instrument_id"] == "en06":
                science_hemisphere_list.append("north")
                target_list[science_frame["request_id"]]["science"][
                    "hemisphere"
                ] = "north"
            elif science_frame["instrument_id"] == "en12":
                science_hemisphere_list.append("south")
                target_list[science_frame["request_id"]]["science"][
                    "hemisphere"
                ] = "south"
            else:
                print("This frame is not generated by FLOYDS.")
            target_list[science_frame["request_id"]]["science"][
                science_frame["OBSTYPE"]
            ] = science_frame["filename"]
            target_list[science_frame["request_id"]]["science"][
                "DAY_OBS"
            ] = science_frame["DAY_OBS"]
            target_list[science_frame["request_id"]]["science"][
                "OBJECT"
            ] = science_frame["OBJECT"]
            # Download the frames
            if (".tar.gz" not in standard_frame["filename"]) or (
                ".fits.fz" not in standard_frame["filename"]
            ):
                [
                    x.append(y)
                    for x, y in zip(
                        [standard_filepath_list, standard_filename_list],
                        download_frame(
                            frame=standard_frame,
                            base_directory=input_folder,
                            no_date=False,
                        ),
                    )
                ]
                if standard_frame["instrument_id"] == "en06":
                    standard_hemisphere_list.append("north")
                    target_list[science_frame["request_id"]]["standard"][
                        "hemisphere"
                    ] = "north"
                elif standard_frame["instrument_id"] == "en12":
                    standard_hemisphere_list.append("south")
                    target_list[science_frame["request_id"]]["standard"][
                        "hemisphere"
                    ] = "south"
                else:
                    print("This frame is not generated by FLOYDS.")
                target_list[science_frame["request_id"]]["standard"][
                    standard_frame["OBSTYPE"]
                ] = standard_frame["filename"]
                target_list[science_frame["request_id"]]["standard"][
                    "DAY_OBS"
                ] = standard_frame["DAY_OBS"]
                target_list[science_frame["request_id"]]["standard"][
                    "OBJECT"
                ] = standard_frame["OBJECT"]

    for x, y in zip(science_filepath_list, science_filename_list):
        print("Science frame: {} is downloaded.".format(x + os.sep + y))

    for x, y in zip(standard_filepath_list, standard_filename_list):
        print("Standard frame: {} is downloaded.".format(x + os.sep + y))

# if working on local files, only works with LCO naming standard for FLOYDS
# north/south is not konwn at this point
# mixed north-south data in the same folder is not allowed
# only accept 3 files of science and 3 files of standard
else:
    # grab all the northern light frames
    filelist_north = glob.glob(
        os.path.join(input_folder_abs_path, "ogg2m001-en06*")
    )
    # grab all the southern light frames
    filelist_south = glob.glob(
        os.path.join(input_folder_abs_path, "coj2m001-en12*")
    )
    print(os.path.join(input_folder_abs_path, "ogg2m001-en06*"))
    if (len(filelist_north) > 0) and (len(filelist_south) > 0):
        raise ValueError("Cannot handle mixed north-south data.")
    elif len(filelist_north) == 6:
        hemisphere = "north"
        filelist = filelist_north
    elif len(filelist_south) == 6:
        hemisphere = "south"
        filelist = filelist_south
    else:
        raise ValueError(
            f"{len(filelist_north)} northern and"
            f" {len(filelist_south)} southern files are provided. Can only"
            " accept 3 science + 3 standard from the same hemisphere."
        )

    for f in filelist:
        ftype = f.split(".")[0].split("-")[-1][0]
        if (
            fits.open(f)[1].header["PROPID"].lower()
            == "FLOYDS standards".lower()
        ):
            if ftype == "e":
                standard_light_frame = f.split(os.path.sep)[-1]
            elif ftype == "w":
                standard_flat_frame = f.split(os.path.sep)[-1]
            elif ftype == "a":
                standard_arc_frame = f.split(os.path.sep)[-1]
            else:
                raise ValueError(
                    f"Unkonwn file type: {ftype}. Should be e, w or a."
                )
        else:
            if ftype == "e":
                science_light_frame = f.split(os.path.sep)[-1]
            elif ftype == "w":
                science_flat_frame = f.split(os.path.sep)[-1]
            elif ftype == "a":
                science_arc_frame = f.split(os.path.sep)[-1]
            else:
                raise ValueError(
                    f"Unkonwn file type: {ftype}. Should be e, w or a."
                )

    target_list = {0: {"science": {}, "standard": {}}}
    target_list[0]["science"]["hemisphere"] = hemisphere
    target_list[0]["science"]["SPECTRUM"] = science_light_frame
    target_list[0]["science"]["LAMPFLAT"] = science_flat_frame
    target_list[0]["science"]["ARC"] = science_arc_frame
    target_list[0]["standard"]["SPECTRUM"] = standard_light_frame
    target_list[0]["standard"]["LAMPFLAT"] = standard_flat_frame
    target_list[0]["standard"]["ARC"] = standard_arc_frame
    print(science_light_frame.split("-")[2])
    target_list[0]["science"]["DAY_OBS"] = science_light_frame.split("-")[2]


yaml_config_list = []

# For each science target
# generate yaml configuration files floyds_target_filename.yaml
for k, v in target_list.items():
    # Star from the template
    with open(yaml_template) as f:
        list_yaml = yaml.load(f)
    # Set the hemisphere north/south
    list_yaml["hemisphere"] = v["science"]["hemisphere"]

    # Modify the names of the frames and output paths
    list_yaml["science_light_frame"] = [v["science"]["SPECTRUM"]]
    list_yaml["science_flat_frame"] = [v["science"]["LAMPFLAT"]]
    list_yaml["science_arc_frame"] = [v["science"]["ARC"]]
    list_yaml["standard_light_frame"] = [v["standard"]["SPECTRUM"]]
    list_yaml["standard_flat_frame"] = [v["standard"]["LAMPFLAT"]]
    list_yaml["standard_arc_frame"] = [v["standard"]["ARC"]]

    list_yaml["target_name"] = target_name
    list_yaml["ra"] = ra
    list_yaml["dec"] = dec
    list_yaml["reducer"] = reducer
    list_yaml["observer"] = observer
    list_yaml["input_folder"] = input_folder
    list_yaml["input_folder"] = input_folder_abs_path
    list_yaml["output_folder"] = output_folder_abs_path
    list_yaml["output_file_name_suffix"] = target_name
    yaml_output_name = "floyds_{}_{}_{}.yaml".format(
        target_name, v["science"]["DAY_OBS"].replace("-", ""), k
    )
    yaml_config_list.append(yaml_output_name)
    if not os.path.exists(output_folder_abs_path):
        os.makedirs(output_folder_abs_path)
    with open(
        os.path.join(output_folder_abs_path, yaml_output_name),
        "w+",
        encoding="ascii",
    ) as f:
        yaml.dump(list_yaml, f)

for yaml_filename in yaml_config_list:
    print("{} will be loaded.".format(yaml_filename))

# run the reduction the reduce_floyds_data.py
for yaml_filename in yaml_config_list:
    os.system(
        f"{sys.executable} {HERE}{os.sep}reduce_floyds_data.py"
        f" {os.path.join(output_folder_abs_path, yaml_output_name)}"
    )
    with open(
        os.path.join(output_folder_abs_path, yaml_output_name), "r"
    ) as stream:
        params = yaml.load(stream)
    outtext = ""
    # Print the data used in a txt file
    outtext += datetime.now().strftime("%Y-%m-%d %H:%M:%S") + os.linesep
    outtext += "Target: {}".format(target_name) + os.linesep
    outtext += "Reduced by: {}".format(reducer) + os.linesep
    outtext += (
        "input_folder_path: {}".format(params["input_folder"]) + os.linesep
    )
    outtext += (
        "Science light frame: {}".format(params["science_light_frame"])
        + os.linesep
    )
    outtext += (
        "Science flat frame: {}".format(params["science_flat_frame"])
        + os.linesep
    )
    outtext += (
        "Science arc frame: {}".format(params["science_arc_frame"])
        + os.linesep
    )
    outtext += (
        "Standard light frame: {}".format(params["standard_light_frame"])
        + os.linesep
    )
    outtext += (
        "Standard flat frame: {}".format(params["standard_flat_frame"])
        + os.linesep
    )
    outtext += (
        "Standard arc frame: {}".format(params["standard_arc_frame"])
        + os.linesep
    )
    outtext += (
        "Intermediate and final data reduction products are saved at: {}".format(
            params["output_folder"]
        )
        + os.linesep
    )
    text_file = open(
        os.path.join(
            output_folder_abs_path, os.path.splitext(yaml_output_name)[0]
        )
        + ".txt",
        "w+",
        encoding="ascii",
    )
    text_file.write(outtext)
    text_file.close()
