from __future__ import division,print_function
from glob import glob
import json
from sys import argv,path
path.append("/home/micha/phd")
try:
    from oscillators import cudamoto_manager as cm
except ImportError:
    import cudamoto_manager as cm

"""
Calculate std and trending of rhist files and write to new file
Usage: arg1 is a globstring of trhist files arg2 is an output directory
"""


globstring = argv[1]

flist = glob(globstring)
out_dir = argv[2]
for idx,fname in enumerate(flist):
    print("File %i of %i (%s)"%(idx,len(flist),fname))
    try:
        cm.reduce_trhist(fname,out_dir)
    except:
        print("unable to process")

