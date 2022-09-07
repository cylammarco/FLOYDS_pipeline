#!/usr/bin/python3
from argparse import ArgumentParser
from collections import OrderedDict
from datetime import timedelta, datetime
import getpass
import json
from multiprocessing.sharedctypes import Value
import os
import sys

import numpy as np
import requests
import ruamel.yaml

from query_lco_archive import get_metadata, download_frame

# Configure the parser
parser = ArgumentParser(
    description="Downloads and Reduce data from archive.lco.global."
)

parser.add_argument(
    "--reducer",
    default="Automated",
    help="Your name to be added to the header.",
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
    "--directory",
    default=None,
    help="Path to store the raw and reduced data products.",
)
parser.add_argument("--login", default=None, help="Path to the login details.")
parser.add_argument(
    "--lco_token",
    default=None,
    help="LCO token. Only used if --login is None. Will ask for one if"
    "needed and nothing is provided.",
)
parser.add_argument(
    "--tns_bot_id",
    default=None,
    help="TNS Bot ID. Only used if --login is None. Will ask for one if"
    "needed and nothing is provided.",
)
parser.add_argument(
    "--tns_bot_name",
    default=None,
    help="TNS Bot name. Only used if --login is None. Will ask for one if"
    "needed and nothing is provided.",
)
parser.add_argument(
    "--tns_token",
    default=None,
    help="TNS token. Only used if --login is None. Will ask for one if"
    "needed and nothing is provided.",
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


# Perform a search on the Transient Name Server
# https://www.wis-tns.org/
def tns_search(search_obj):
    search_url = url_tns_api + "/search"
    headers = {
        "User-Agent": 'tns_marker{"tns_id": "'
        + str(TNS_BOT_ID)
        + '", "type": "bot", "name": "'
        + TNS_BOT_NAME
        + '"}'
    }
    json_file = OrderedDict(search_obj)
    search_data = {"api_key": TNS_API_KEY, "data": json.dumps(json_file)}
    response = requests.post(search_url, headers=headers, data=search_data)
    return response


# Download the info of the object of interest from the TNS
def get(get_obj):
    get_url = url_tns_api + "/object"
    headers = {
        "User-Agent": 'tns_marker{"tns_id": "'
        + str(TNS_BOT_ID)
        + '", "type": "bot", "name": "'
        + TNS_BOT_NAME
        + '"}'
    }
    json_file = OrderedDict(get_obj)
    get_data = {"api_key": TNS_API_KEY, "data": json.dumps(json_file)}
    response = requests.post(get_url, headers=headers, data=get_data)
    return response


# Load the login information: TNS & SNEx
if args.login is not None:

    with open(args.login) as f:
        login_yaml = yaml.load(f)

# choose directory to store the output
base_directory = args.directory
target_name = args.target_name
reducer = args.reducer

if (base_directory is None) & (target_name is None):

    raise ValueError("Either base_directory or target_name has to be provided")

elif base_directory is None:

    base_directory = target_name

elif target_name is None:

    if (args.ra is None) or (args.dec is None):

        raise ValueError(
            "Either target_name or (ra & dec) has to be provided."
        )

    target_name = base_directory
    ra = args.ra
    dec = args.dec

else:

    print("target name is provided, input ra and dec are ignored.")

    # get the target position if only a name is provided
    TNS = "www.wis-tns.org"
    url_tns_api = "https://" + TNS + "/api/get"

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
    get_obj_response = get(get_obj)

    ra = get_obj_response.json()["data"]["reply"]["radeg"]
    dec = get_obj_response.json()["data"]["reply"]["decdeg"]


# get the token
# LCO token: needed for private data
# TNS token: needed for name resolving to (ra, dec)
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
                    datetime.fromisoformat(meta["observation_date"][:-1])
                )
            _time_diff = [i - _date for i in _time]
            min_time_idx = np.argmin(np.abs(np.array(_time_diff)))
        elif len(_standard_metadata_day_obs) == 1:
            min_time_idx = 0
            _time.append(datetime.fromisoformat(meta["observation_date"][:-1]))
        else:
            # Get all the timestamps
            for meta in _standard_metadata:
                print(meta["DAY_OBS"])
                _time.append(
                    datetime.fromisoformat(meta["observation_date"][:-1])
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
            start=(_time[min_time_idx] - timedelta(minutes=30)).isoformat(),
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
for i in range(len(science_metadata)):
    science_frame = science_metadata[i]
    standard_frame = standard_metadata[i]
    if (".tar.gz" not in science_frame["filename"]) & (
        science_frame["request_id"] in request_id_science
    ):
        [
            x.append(y)
            for x, y in zip(
                [science_filepath_list, science_filename_list],
                download_frame(
                    frame=science_frame,
                    base_directory=base_directory,
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
        if ".tar.gz" not in standard_frame["filename"]:
            [
                x.append(y)
                for x, y in zip(
                    [standard_filepath_list, standard_filename_list],
                    download_frame(
                        frame=standard_frame,
                        base_directory=base_directory,
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
    list_yaml["science_light_frame"] = [
        os.path.join(
            base_directory,
            v["science"]["SPECTRUM"].split("-")[2],
            v["science"]["SPECTRUM"],
        )
    ]
    list_yaml["science_flat_frame"] = [
        os.path.join(
            base_directory,
            v["science"]["LAMPFLAT"].split("-")[2],
            v["science"]["LAMPFLAT"],
        )
    ]
    list_yaml["science_arc_frame"] = [
        os.path.join(
            base_directory,
            v["science"]["ARC"].split("-")[2],
            v["science"]["ARC"],
        )
    ]
    list_yaml["standard_light_frame"] = [
        os.path.join(
            base_directory,
            v["standard"]["SPECTRUM"].split("-")[2],
            v["standard"]["SPECTRUM"],
        )
    ]
    list_yaml["standard_flat_frame"] = [
        os.path.join(
            base_directory,
            v["standard"]["LAMPFLAT"].split("-")[2],
            v["standard"]["LAMPFLAT"],
        )
    ]
    list_yaml["standard_arc_frame"] = [
        os.path.join(
            base_directory,
            v["standard"]["ARC"].split("-")[2],
            v["standard"]["ARC"],
        )
    ]
    list_yaml["target_name"] = target_name
    list_yaml["reducer"] = reducer
    list_yaml["output_folder_path"] = os.path.join(
        base_directory, v["science"]["DAY_OBS"].replace("-", ""), "output"
    )
    list_yaml["output_file_name_suffix"] = base_directory
    yaml_output_name = "floyds_{}_{}_{}.yaml".format(
        target_name, v["science"]["DAY_OBS"].replace("-", ""), k
    )
    yaml_config_list.append(yaml_output_name)
    with open(
        yaml_output_name,
        "w+",
    ) as f:
        yaml.dump(list_yaml, f)

for yaml_filename in yaml_config_list:
    print("{} will be loaded.".format(yaml_filename))

# run the reduction the reduce_floyds_data.py
for yaml_filename in yaml_config_list:

    os.system(
        "{} {}{}reduce_floyds_data.py {}".format(
            sys.executable, HERE, os.sep, yaml_filename
        )
    )
    with open(yaml_filename, "r") as stream:
        params = yaml.load(stream)
    outtext = ""
    # Print the data used in a txt file
    outtext += datetime.now().strftime("%Y-%m-%d %H:%M:%S") + os.linesep
    outtext += "Target: {}".format(target_name) + os.linesep
    outtext += "Reduced by: {}".format(reducer) + os.linesep
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
            params["output_folder_path"]
        )
        + os.linesep
    )
    text_file = open(
        os.path.join(
            base_directory,
            "_".join(params["output_folder_path"].split(os.path.sep)[:-1])
            + ".txt",
        ),
        "w+",
    )
    text_file.write(outtext)
    text_file.close()
