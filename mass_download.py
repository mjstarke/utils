import shutil
from urllib.request import urlopen
from contextlib import closing

"""
1. Save the list of file names, one per line, to some local file.
2. Set LISTFILE to point to that file.
3. If the file names do not include the full path, set BASELINK to the part that is omitted.
4. Set EXCERPT_START to something at the start of the file name, and EXCERPT_END to something at the end.
  This allows the list of files to be copied from a table in whole, even if additional columns are also copied.
5. Set DOWNLOAD_TARGET to hwere you want to download files locally.
"""

LISTFILE = "mass_download.txt"
BASELINK = 'http://noaa-nexrad-level2.s3.amazonaws.com/2018/06/29/KDLH/'
EXCERPT_START = "K"
EXCERPT_END = "V06"
DOWNLOAD_TARGET = "download"


def excerpt(s, start, end):
    return s[s.index(start):s.index(end) + len(end)]


targets = []
with open(LISTFILE, "r") as f:
    for line in f:
        if (EXCERPT_START in line) and (EXCERPT_END in line):
            targets.append(BASELINK + excerpt(line, EXCERPT_START, EXCERPT_END))

ex = []

for t in range(len(targets)):
    try:
        target = targets[t]
        print()
        "({:3} of {:3}) {:28}  {}".format(t + 1, len(targets), "Starting download:", target)
        with closing(urlopen(target)) as r:
            with open(DOWNLOAD_TARGET + "/" + target.split("/")[-1], 'wb') as f:
                shutil.copyfileobj(r, f)
        print()
        "({:3} of {:3}) {:28}".format(t + 1, len(targets), "Successfully downloaded.")
    except Exception as e:
        ex.append(e)
        print()
        "({:3} of {:3}) {:28}".format(t + 1, len(targets), "Download failed:")
        print()
