#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  : Nishant Kumar
# Created Date: 25/02/2022
# ---------------------------------------------------------------------------
""" 
Convert POD modes to VTK
"""

import numpy as np
import pandas as po
import sys
from tqdm import tqdm
from pathlib import Path
from pyevtk.hl import pointsToVTK

# import time as clock
# start_time = clock.time()

# ---------------------------------------------------------------------------
# UTILITY FUNCTION(S)
# ---------------------------------------------------------------------------
#- Count number of lines in a file
#- See: https://stackoverflow.com/a/845081/7473705
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#- Check and make directory path
def make_dir(dirpath):
    Path(dirpath).mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# SUBFUNCTION(S)
# ---------------------------------------------------------------------------
class convertModeToVTK():
    def __init__(self, pointFile, modeDir, nModes):
        self.pointFile = pointFile
        self.modeDir = modeDir
        self.nModes = nModes
        self.varSize = 3

    def readData(self):
        #- Coordinates
        if self.pointFile.endswith('.dat'):
            useInd = [1,2,3]
        elif self.pointFile.endswith('.xy'):
            useInd = [0,1,2]
        pts = po.read_csv(self.pointFile, 
                          delim_whitespace=True, 
                          header=None, 
                          names=['x', 'y', 'z'], 
                          usecols=useInd)
        self.x = pts[['x']].to_numpy().squeeze()
        self.y = pts[['y']].to_numpy().squeeze()
        self.z = pts[['z']].to_numpy().squeeze()

        #- Size
        self.MM = len(self.x)

        self.x = self.x.reshape(self.MM, 1, order='C')
        self.y = self.y.reshape(self.MM, 1, order='C')
        self.z = self.z.reshape(self.MM, 1, order='C')

        #- Read the reconstructed field
        print('\nReading modes from binary files...\n')
        self.mode = np.fromfile(self.modeDir+'/mode.bin', dtype=float)
        self.mode = self.mode.reshape(self.nModes, self.varSize*self.MM)

    def saveVTK(self):
        vtkDir = self.modeDir+"/VTK/"
        make_dir(vtkDir)
        print('\nWriting the modes in VTK files...\n')
        MM = self.MM
        for i in tqdm(range(self.nModes)):
            modei = (self.mode[i, :MM], self.mode[i, MM:2*MM], self.mode[i, 2*MM:3*MM])
            pointsToVTK(vtkDir+"mode_%s"%(str(i)), 
                    self.x, self.y, self.z, data={"mode": modei})

# ---------------------------------------------------------------------------
# MAIN FUNCTION
# ---------------------------------------------------------------------------
def main():
    #- Read arguments
    #- NOTE: sys.argv[0] is always the script name
    pointFile = sys.argv[1]
    modeDir = sys.argv[2]
    nModes = int(sys.argv[3])

    p = convertModeToVTK(pointFile, modeDir, nModes)
    p.readData()
    p.saveVTK()

# ---------------------------------------------------------------------------
# COMMAND LINE EXECUTION
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    main()