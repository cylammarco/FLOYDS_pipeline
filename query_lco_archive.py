import os
import requests
import re

# This script is based on
# https://github.com/LCOGT/lcogtsnpipe/blob/master/trunk/bin/LCOGTingest.py

# Refer also to this page for more information on the parameters/arguments:
# https://lco.global/documentation/archive-documentation/


def get_metadata(authtoken={}, limit=None, logger=None, **kwargs):
    """
    Get the list of files meeting criteria in kwargs

    """

    url = "https://archive-api.lco.global/frames/?" + "&".join(
        [
            key + "=" + str(val)
            for key, val in kwargs.items()
            if val is not None
        ]
    )
    url = url.replace("False", "false")
    url = url.replace("True", "true")
    if logger is not None:
        logger.info(url)
    else:
        print(url)

    response = requests.get(url, headers=authtoken, stream=True).json()
    frames = response["results"]
    while response["next"] and (limit is None or len(frames) < limit):
        if logger is not None:
            logger.info(response["next"])
        else:
            print(response["next"])
        response = requests.get(
            response["next"], headers=authtoken, stream=True
        ).json()
        frames += response["results"]
    return frames[:limit]


def download_frame(frame, base_directory, no_date=False):
    """
    Download a single image from the LCOGT archive

    Parameters
    ==========
    frame: dict
        Dictionary containing the response from requests.
    base_directory: str
        The relative/full path of where the files are to be stored.
    no_date: boolean
        Set to True to store everything in base_directory, else, the files
        will be stored in folders named by the date of the night of the
        observation.

    """

    filename = frame["filename"]
    dayobs = re.search("(20\d\d)(0\d|1[0-2])([0-2]\d|3[01])", filename).group()
    if no_date:
        filepath = os.path.join(base_directory)
    else:
        filepath = os.path.join(base_directory, dayobs)

    if not os.path.isdir(filepath):
        os.makedirs(filepath)

    filename = frame["filename"]
    with open(os.path.join(filepath, filename), "wb") as f:
        f.write(requests.get(frame["url"]).content)

    return filepath, filename
