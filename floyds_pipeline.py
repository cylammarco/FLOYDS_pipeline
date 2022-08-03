#!/usr/bin/python3
from query_lco_archive import get_metadata, download_frame
import os
from collections import OrderedDict
import requests
import json
from datetime import date, timedelta, datetime
import numpy as np
import ruamel.yaml
import sys

yaml = ruamel.yaml.YAML()
yaml.preserve_quotes = True


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


base_directory = "AT2022jdf"
target_name = "2022jdf"
# choose directory to store the output


# get the token
# LCO token: needed for private data
# TNS token: needed for name resolving to (ra, dec)
lco_token = ""

if lco_token is None:
    authtoken = {}
else:
    authtoken = {"Authorization": "Token {}".format(lco_token)}


# get the target position if only a name is provided
TNS = "www.wis-tns.org"
url_tns_api = "https://" + TNS + "/api/get"

TNS_BOT_ID = ""
TNS_BOT_NAME = ""
TNS_API_KEY = ""

search_obj = [("objname", target_name)]
tns_search_response = tns_search(search_obj)

# assuming only 1 result...
get_obj = [("objid", tns_search_response.json()["data"]["reply"][0]["objid"])]
get_obj_response = get(get_obj)

ra = get_obj_response.json()["data"]["reply"]["radeg"]
dec = get_obj_response.json()["data"]["reply"]["decdeg"]


# get the filelist
# en06 is FLOYDS North at the Haleakala Observatory (OGG)
# en12 is FLOYDS South at the Siding Spring Observatory (COJ)
# the search has to be done seprately

start_mjd = "1900-01-01"
end_mjd = "2100-12-31"

science_metadata = []
standard_metadata = []
request_id = []
obstype = []
day_obs = []

# first north (en06) then south (en12)
for instrume in ["en06", "en12"]:

    # Get the metadata
    _science_metadata = get_metadata(
        authtoken=authtoken,
        limit=1000,
        INSTRUME=instrume,
        start=start_mjd,
        end=end_mjd,
        RLEVEL=0,
        covers="POINT({} {})".format(ra, dec),
    )

    for metadata in _science_metadata:
        request_id.append(metadata["request_id"])
        obstype.append(metadata["OBSTYPE"])
        day_obs.append(metadata["DAY_OBS"])

    # Get the request_id that contains useful spectral data
    request_id_science = [
        i for i, j in zip(request_id, obstype) if j == "SPECTRUM"
    ]
    request_id_standard = []

    science_metadata += [
        v for v in _science_metadata if v["request_id"] in request_id_science
    ]

    for day in [i for i, j in zip(day_obs, obstype) if j == "SPECTRUM"]:

        _day = datetime.fromisoformat(day)
        day_range = 0

        _standard_metadata = []
        while _standard_metadata == []:
            # Get the standard star
            _standard_metadata = get_metadata(
                authtoken=authtoken,
                limit=1000,
                INSTRUME=instrume,
                start=(_day - timedelta(days=day_range)).isoformat(),
                end=(_day + timedelta(days=day_range)).isoformat(),
                PROPID="FLOYDS standards",
                OBSTYPE="SPECTRUM",
                RLEVEL=0,
            )
            day_range += 2

        # If there are more than 1 frame returned, keep the standard with the
        # least time difference
        _time = []
        if len(_standard_metadata) > 1:
            # Get all the timestamps
            for meta in _standard_metadata:
                _time.append(
                    datetime.fromisoformat(meta["observation_date"][:-1])
                )

            _time_diff = [i - _day for i in _time]
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
        day_range = 0
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
        day_range = 0
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
        if science_frame["instrument_id"] == "en06":
            science_hemisphere_list.append("north")
        elif science_frame["instrument_id"] == "en12":
            science_hemisphere_list.append("south")
        else:
            print("This frame is not generated by FLOYDS.")

        # Arrange into the target list
        if science_frame["request_id"] not in target_list.keys():
            target_list[science_frame["request_id"]] = {
                "science": {},
                "standard": {},
            }
        target_list[science_frame["request_id"]]["science"][
            science_frame["OBSTYPE"]
        ] = science_frame["filename"]
        target_list[science_frame["request_id"]]["science"][
            "DAY_OBS"
        ] = science_frame["DAY_OBS"]
        target_list[science_frame["request_id"]]["science"][
            "OBJECT"
        ] = science_frame["OBJECT"]

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
            elif standard_frame["instrument_id"] == "en12":
                standard_hemisphere_list.append("south")
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
    with open("floyds_template.yaml") as f:
        list_yaml = yaml.load(f)
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


# use backup arc and flat (north-south dependent) if any of the science or standard calibration frames are missing


# run the reduction the reduce_floyds_data.py
for yaml_filename in yaml_config_list:
    os.system(
        "{} reduce_floyds_data.py {}".format(sys.executable, yaml_filename)
    )
