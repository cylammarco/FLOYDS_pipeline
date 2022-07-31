from argparse import ArgumentParser
import logging
import os
import requests
import re

# This script is based on
# https://github.com/LCOGT/lcogtsnpipe/blob/master/trunk/bin/LCOGTingest.py

# Refer also to this page for more information on the parameters/arguments:
# https://lco.global/documentation/archive-documentation/


def get_metadata(authtoken={}, limit=None, **kwargs):
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
    logger.info(url)

    response = requests.get(url, headers=authtoken, stream=True).json()
    frames = response["results"]
    while response["next"] and (limit is None or len(frames) < limit):
        logger.info(response["next"])
        response = requests.get(
            response["next"], headers=authtoken, stream=True
        ).json()
        frames += response["results"]
    return frames[:limit]


def download_frame(frame, base_directory, no_date):
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


logger = logging.getLogger()

parser = ArgumentParser(description="Downloads data from archive.lco.global")

# True / False
parser.add_argument("-D", "--no-date", action="store_true")
parser.add_argument(
    "--orac",
    action="store_true",
    help="get products from the old ORAC-DR pipeline",
)
parser.add_argument(
    "--public", action="store_true", help="include public data"
)


# choices
parser.add_argument(
    "-S",
    "--site",
    choices=["bpl", "coj", "cpt", "elp", "lsc", "ogg", "sqa", "tfn"],
    default=None,
)
parser.add_argument(
    "-T",
    "--telescope",
    choices=["0m4a", "0m4b", "0m4c", "0m8a", "1m0a", "2m0a"],
    default=None,
)
parser.add_argument(
    "-f",
    "--filter",
    choices=["up", "gp", "rp", "ip", "zs", "U", "B", "V", "R", "I"],
    default=None,
)
parser.add_argument(
    "-t",
    "--obstype",
    choices=[
        "ARC",
        "BIAS",
        "CATALOG",
        "DARK",
        "EXPERIMENTAL",
        "EXPOSE",
        "LAMPFLAT",
        "SKYFLAT",
        "SPECTRUM",
        "STANDARD",
    ],
    default=None,
)
parser.add_argument(
    "-r", "--reduction", choices=["raw", "quicklook", "reduced"], default=None
)

# free form
parser.add_argument("-d", "--directory", default=None)
parser.add_argument("-a", "--token", default=None)
parser.add_argument(
    "-l",
    "--limit",
    type=int,
    help="maximum number of frames to return",
    default=10000,
)
parser.add_argument("-I", "--instrument", default=None)
parser.add_argument(
    "-P", "--proposal", help="proposal ID (PROPID in the header)", default=None
)
parser.add_argument("-n", "--name", help="target name", default=None)
parser.add_argument("-s", "--start", help="start date", default=None)
parser.add_argument("-e", "--end", help="end date", default=None)
parser.add_argument(
    "-c",
    "--coords",
    nargs=2,
    help="target coordinates in degrees, space separated",
    default=None,
)


args = parser.parse_args()

if args.reduction == "raw":
    rlevel = 0
elif args.reduction == "quicklook" and args.orac:
    rlevel = 10
elif args.reduction == "quicklook":
    rlevel = 11
elif args.reduction == "reduced" and args.orac:
    rlevel = 90
elif args.reduction == "reduced":
    rlevel = 91
else:
    rlevel = None


if args.directory is None:

    directory = args.name

else:

    directory = args.directory

# Get your token at https://observe.lco.global/accounts/profile
if args.token is None:
    authtoken = {}
else:
    authtoken = {"Authorization": "Token {}".format(args.token)}

# 0718806305b4e990d740c782d77762665e02dadc


frames = get_metadata(
    authtoken,
    limit=args.limit,
    SITEID=args.site,
    TELID=args.telescope,
    INSTRUME=args.instrument,
    FILTER=args.filter,
    PROPID=args.proposal,
    OBJECT=args.name,
    start=args.start,
    end=args.end,
    OBSTYPE=args.obstype,
    RLEVEL=rlevel,
    public=args.public,
    covers="POINT({} {})".format(*args.coords) if args.coords else None,
)


for frame in frames:
    if ".tar.gz" not in frame["filename"]:
        filepath, filename = download_frame(
            frame=frame, base_directory=directory, no_date=args.no_date
        )
