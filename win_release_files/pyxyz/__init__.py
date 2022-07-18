import os, sys


# if "LD_LIBRARY_PATH" not in os.environ or mypath not in os.environ["LD_LIBRARY_PATH"]:
    # print("[WARNING] In case you don't have gslcblas installed and get ImportError, call 'export LD_LIBRARY_PATH={}:$LD_LIBRARY_PATH'".format(mypath))

mypath = os.path.dirname(os.path.realpath(__file__))
if mypath not in sys.path:
    sys.path.insert(0, mypath)
os.add_dll_directory(mypath)

from confpool import Confpool
