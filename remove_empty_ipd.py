"""
Given a folder of IPD files, remove the IPD files with only headers and no data.
"""

import os
import sys





def remove_empty_ipd(input_folder):
    for file in os.listdir(input_folder):
        if file.endswith(".ipd1.csv"):
            with open(os.path.join(input_folder, file), "r") as f:
                lines = f.readlines()
                if len(lines) < 10:
                    print(f"Removing {file} with {len(lines)} lines")
                    os.remove(os.path.join(input_folder, file))


input_folder = os.path.join(sys.argv[1], 'ipd')
remove_empty_ipd(input_folder)