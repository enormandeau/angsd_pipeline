#!/usr/bin/env python3
"""Filter SNPS from the output of 03_saf_maf_gl_depth_all.sh

Usage:
    <program> angsd_output_folder filtered_folder window_size

The script will read all needed files in the angsd_output_folder and use the
information in all of them to filter the SNPs in one pass.

Files found in the angsd_output_folder:

*_maf0.01_pctind0.5_maxdepth20.arg
*_maf0.01_pctind0.5_maxdepth20.beagle.gz    # used
*_maf0.01_pctind0.5_maxdepth20.counts.gz    # used
*_maf0.01_pctind0.5_maxdepth20.depthGlobal
*_maf0.01_pctind0.5_maxdepth20.depthSample
*_maf0.01_pctind0.5_maxdepth20.mafs         # used
*_maf0.01_pctind0.5_maxdepth20.pos.gz
*_maf0.01_pctind0.5_maxdepth20.saf.gz
*_maf0.01_pctind0.5_maxdepth20.saf.idx
*_maf0.01_pctind0.5_maxdepth20.saf.pos.gz
"""

# Modules
from collections import defaultdict
import gzip 
import sys
import os

# Function
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

# Parse user input
try:
    angsd_output_folder = sys.argv[1]
    filtered_folder = sys.argv[2]
    window_size = int(sys.argv[3])
except:
    print(__doc__)
    sys.exit(1)

# Open input file handles
files = os.listdir(angsd_output_folder)
input_files = {}

input_files["beagle"] = myopen(
        os.path.join(
            angsd_output_folder,
            [f for f in files if f.endswith("beagle.gz")][0]
            )
        )

input_files["counts"] = myopen(
        os.path.join(
            angsd_output_folder,
            [f for f in files if f.endswith("counts.gz")][0]
            )
        )

input_files["mafs"] = myopen(
        os.path.join(
            angsd_output_folder,
            [f for f in files if f.endswith("mafs")][0]
            )
        )

# Create output folder
if not os.path.exists(filtered_folder):
    os.mkdir(filtered_folder)

# Open output file handles
files = os.listdir(angsd_output_folder)
output_files = {}

output_files["beagle"] = myopen(
        os.path.join(
            filtered_folder,
            [f for f in files if f.endswith("beagle.gz")][0]
            ),
        "wt"
        )

output_files["counts"] = myopen(
        os.path.join(
            filtered_folder,
            [f for f in files if f.endswith("counts.gz")][0]
            ),
        "wt"
        )

output_files["mafs"] = myopen(
        os.path.join(
            filtered_folder,
            [f for f in files if f.endswith("mafs")][0]
            ),
        "wt"
        )

# SNP containers
left_snps = []
central_snp = None
rigtht_snps = []
surrounding = []

# Write headers to output files
for f in input_files:
    output_files[f].write(input_files[f].readline())

# Read lines from all 3 files one at a time
for mafs_line in input_files["mafs"]:

    # Position info (scaffold<str>, position<int>)
    pos = mafs_line.strip().split()[:2]
    pos = (pos[0], int(pos[1]))

    # Keep infos about positions in 'queue' in order to keep as little
    # information as needed to treat a given SNP. For example, we need only the
    # 30 positions on each side of a SNP to decide if we are keeping it.

    # Format mafs for later use
    mafs = mafs_line.strip().split()
    mafs = [
            mafs[0],
            int(mafs[1]),
            mafs[2],
            mafs[3],
            mafs[4],
            float(mafs[5]),
            int(mafs[6])
            ]

    # Format beagle for later use
    beagle = input_files["beagle"].readline().strip().split()

    # Format count for later use
    counts = input_files["counts"].readline().strip().split()
    counts = [int(x) for x in counts]

    info = (pos, mafs, beagle, counts)

    if not central_snp:
        central_snp = info

    elif len(rigtht_snps) < window_size:
        rigtht_snps.append(info)
        print(len(rigtht_snps))

    else:
        pass
        # - Filter
        # - move SNPs from right to left

# Close file handles
for f in input_files:
    input_files[f].close()

for f in output_files:
    output_files[f].close()
