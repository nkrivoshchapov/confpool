import os, sys

mypath = os.path.dirname(os.path.realpath(__file__))
if "LD_LIBRARY_PATH" not in os.environ or mypath not in os.environ["LD_LIBRARY_PATH"]:
    print("[WARNING] In case you don't have gslcblas installed, call 'export LD_LIBRARY_PATH={}:$LD_LIBRARY_PATH'".format(mypath))

if mypath not in sys.path:
    sys.path.insert(0, mypath)

from confpool import *
