import time
import os
from Bio import SeqIO
import multiprocessing as mp
import numpy as np
import platform

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
SYSTEM_NM = platform.system()

if SYSTEM_NM == 'Linux':
    # REAL
    REF_DIR = "../hg38/"
else:
    # DEV
    WORK_DIR = "D:/000_WORK/Ramu_JangHyeWon/20201103/WORK_DIR/"
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"


TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
############### end setting env #################

def test():
    util = Util.Utils()
    pass

if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    test()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))