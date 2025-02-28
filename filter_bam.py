import pysam
import numpy as np
import os


from pbcore.io import AlignmentSet
import numpy as np
import os
import pandas as pd
from collections import defaultdict
# from scipy.stats import norm
import time
import sys
import logging

# Input and output BAM files
input_bam = "/home/shuaiw/methylation/data/borg/all_borg/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam"  # Replace with your BAM file
output_bam = "/home/shuaiw/methylation/data/borg/all_borg/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.filter.bam"
input_bam = "/home/shuaiw/methylation/data/borg/break_contigs2/XRSBK_20221007_S64018_PL100268287-1_C01.align.ccs.bam"
aln = AlignmentSet(input_bam)
# Define the int32 maximum value
INT32_MAX = np.iinfo(np.int32).max  # 2147483647

# # Open BAM file for reading and writing
i = 0
with pysam.AlignmentFile(input_bam, "rb") as bam_in:
    header = bam_in.header.to_dict()
    # Remove RG entries with the specified IDs
    ## replace the RG with the specified ID

    for rg in header['RG']:
        print (rg.get('ID'))
    print (header['RG'])

# with pysam.AlignmentFile(input_bam, "rb") as bam_in, \
#      pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:

#     for record in bam_in:
#         ## skip if the read is not mapped
#         if record.is_unmapped:
#             continue
#         # Extract Read Group (RG) tag
#         rg_tag = record.get_tag("RG") if record.has_tag("RG") else None
#         # print (f"Read {record.query_name} has RG {rg_tag}")

#         ## check if the rg_tag is 'cb4d472d-100C60F6' or 'cb4d472d'
#         if rg_tag == 'cb4d472d-100C60F6' or rg_tag == 'cb4d472d':
#             print (f"Read {record.query_name} has RG {rg_tag}, skipping")
#             continue

#         if rg_tag:
#             try:
#                 rg_int = int(rg_tag.split("/")[0], 16)  # Convert from hex to int
#                 if rg_int > INT32_MAX:
#                     print (f"Read {record.query_name} has RG {rg_int}, skipping, {rg_tag}")
#                     continue  # Skip this read if RG is too large
#                 rg_int = np.int32(rg_int)  # Convert to int32
#             except ValueError:
#                 pass  # Skip if RG is not a valid hex number

#         # Write the valid read to the output BAM
#         bam_out.write(record)
#         i += 1
        # if i > 1:
        #     break

# # print(f"Filtered BAM file saved as: {output_bam}")
# os.system(f"samtools index {output_bam}")  # Index the output BAM file
# os.system(f"/home/shuaiw/smrtlink/pbindex {output_bam}")  # Index the output BAM file

# ## read the output bam file using pbcore
aln = AlignmentSet(output_bam)

