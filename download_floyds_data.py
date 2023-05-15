from argparse import ArgumentParser
import logging
import os

from query_lco_archive import get_metadata, download_frame

# This script is based on
# https://github.com/LCOGT/lcogtsnpipe/blob/master/trunk/bin/LCOGTingest.py

# Refer also to this page for more information on the parameters/arguments:
# https://lco.global/documentation/archive-documentation/

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

filepath_list = []
filename_list = []
for frame in frames:
    if ".tar.gz" not in frame["filename"]:
        [
            x.append(y)
            for x, y in zip(
                [filepath_list, filename_list],
                download_frame(
                    frame=frame, base_directory=directory, no_date=args.no_date
                ),
            )
        ]

for x, y in zip(filepath_list, filename_list):
    print("{} is downloaded.{}".format(x + os.sep + y, os.linesep))
