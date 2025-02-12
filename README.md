# Methy

https://github.com/PacificBiosciences/pbcore/issues/127
```
  File "/home/shuaiw/miniconda3/envs/methy/lib/python3.12/site-packages/pbcore/io/align/BamIO.py", line 85, in _loadReadGroupInfo
    rgID = rgAsInt(rg["ID"])
           ^^^^^^^^^^^^^^^^^
  File "/home/shuaiw/miniconda3/envs/methy/lib/python3.12/site-packages/pbcore/io/align/_BamSupport.py", line 66, in rgAsInt
    return np.int32(int(rgIdString.split("/")[0], 16))
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
```
OverflowError: Python integer 3367457666 out of bounds for int32