#!/usr/bin/env python2.7
#
#  plot_edf_wavelets.py
#
#  Created by William Bosl on 2/15/2017.
#  Copyright (c) 2017 MindLight Medical, Inc. All rights reserved.
#


import numpy as np
import pyedflib
import sys

#--------------------------------------------------------
# Usage instructions
#--------------------------------------------------------
def get_commandline(argv):
    argc = len(argv)
    if argc < 2:
        print ("Usage: python read_edf.py anyfile.edf")
        exit()
    else:
        filename = argv[1]
        dicom_filename = filename + ".dicom"
    return filename, dicom_filename

#--------------------------------------------------------
# Read in the data and pull out the data
#--------------------------------------------------------
def read_edf(filename):
    f = pyedflib.EdfReader(filename)
    
    # Number of channels in this file
    nch = f.signals_in_file
    
    # Names of each channel in this file
    channelNames = f.getSignalLabels()

    # Sampling rate
    srate = f.getSampleFrequency(0)
    
    # These are the time series associated with each channel name
    data = np.zeros((nch, f.getNSamples()[0]))
    for i in np.arange(nch):
        data[i, :] = f.readSignal(i)
        
    # Number of points in each time series
    n = data.shape[1]

    # Print some information 
    print ("data size: ", data.shape)
    print ("channels = ", channelNames)
    print ("srate =     ", srate)
    print ("Number of seconds in this file = ", n/srate )
    
    return(nch, channelNames, srate, data)

#--------------------------------------------------------
# Read in the data and pull out the data
#--------------------------------------------------------
def write_DICOM_eeg(dicom_filename):
    print("\n")
    print("Writing data to DICOM file: ", dicom_filename)

#--------------------------------------------------------
# main program - called when this file is executed as
# a stand-alone program
#--------------------------------------------------------
if __name__ == "__main__":    
    
    # Read the output filename if given
    filename, dicom_filename = get_commandline(sys.argv)
    
    # Extract information in the EDF file
    (nch, channelNames, srate, data) = read_edf(filename) 
    
    # Now do something with the EDF file! 
    write_DICOM_eeg(dicom_filename)

