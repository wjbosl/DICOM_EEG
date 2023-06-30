"""read_edfD.py:

__author__ = "William Bosl"
__copyright__ = "Copyright 2022, William J. Bosl"
__credits__ = ["William Bosl"]
__license__ = All rights reserved by William J. Bosl
__version__ = "1.0.0"
__maintainer__ = "William Bosl"
__email__ = "william.bosl@childrens.harvard.edu"
__status__ = "Initial submission"

This file opens a discontinous EDF+D file and:
    - identifies the time stamps for every continous segment
    - writes out all of the header information and segment information
    - optionally writes new EDF+C files with continuous segments
"""

import sys
import numpy as np
import datetime

def filter_nonprintable(text):
    import itertools
    # Use characters of control category
    nonprintable = itertools.chain(range(0x00,0x20),range(0x7f,0xa0))
    # Use translate to remove all non-printable characters
    return text.translate({character:None for character in nonprintable})

# Define the EDF file structure
def file_spec(f):
	size = f.count('size', 'text', 2)  # A two-byte unsigned integer
	f.bytes('text', size)  # A variable number of bytes


#-----------------------------------------
#
#-----------------------------------------
def parse_commandline(argv):

    # Set defaults, then read inputs
    rw = "read"
    argc = len(argv)
    missing_i = True
    if "-i" in argv: missing_i = False
    if argc < 3 or missing_i:
        print("Usage: > python read_edfD.py -i filename.edf <-rw write or read>")
        print("Example: > python read_edfD.py -i somefile.edf -rw read")
        exit()
    else:
        for i in range(1,argc,2):
            option = argv[i]
            value = argv[i+1]
            if option == '-i':  # input filename
                filename = value
            elif option == '-rw':
                rw = value

    return filename, rw


#-----------------------------------------
# Add seconds to date/time
#-----------------------------------------
def add_seconds(date, time, secs):
    # Convert the strings to a datetime object
    a = range(3)
    ymd = date.split('.')
    y,m,d = [int(ymd[i]) for i in a]
    hms = time.split('.')
    h,min,s = [int(hms[i]) for i in a]
    dt = datetime.datetime(y,m,d,h,min,s)
    
    print("old date, time = ", date, time,"; dt = ", secs)


    # Add the given number of seconds
    new_dt = dt + datetime.timedelta(seconds=secs)
    
    # Convert the datetime object back to strings
    new_date = new_dt.date().strftime("%2Y.%m.%d")
    new_time = new_dt.time().strftime("%H.%M.%S")
    
    return new_date, new_time

#-----------------------------------------
# Read the edf file and print some information
#-----------------------------------------
def read_edfD(filename):
    
    print("reading data from ", filename)
    
    # Create a record to hold the size of the header and data 
    # records and location of continuous segments
    record = {}
    record["nhead"] = 0 
    record["nrecords"] = []
    record["nr"] = 0
    record["ns"] = 0
    record["start"] = [] 
    record["nheadbytes"] = 0

    # Read the file and print the text field
    file = open(filename, "rb")
    
    # Get header information
    version = file.read(8).decode('ascii')
    patient = file.read(80).decode('ascii')
    localrecID = file.read(80).decode('ascii')
    startdate = file.read(8).decode('ascii')
    starttime = file.read(8).decode('ascii')
    record["start"].append([startdate, starttime])
    print("Patient field: ", patient)
    print("startdate: ", startdate)
    print("starttime: ", starttime)

    nheadbytes = file.read(8).decode('ascii')
    reserved1 = file.read(44).decode('ascii')
    print("reserved1 = ", reserved1)
    nr = int(file.read(8).decode('ascii'))  # number of records
    rec_in_secs = float(file.read(8).decode('ascii'))
    ns = int(file.read(4).decode('ascii'))  # number of signals
    nhead = int(nheadbytes)
    record["nhead"] = nhead
    record["nr"] = nr
    record["ns"] = ns
    
    edf_type = reserved1
    print("EDF+ type: ", edf_type)
    print("version: ", version)
    print("nr, ns, rec_in_secs: ", nr, ns, rec_in_secs)

    #print("Step 1, position = ", file.tell())
   
    # Read header data for each signal
    new_ch_names = []
    transducer = []
    physical_dim = []
    physical_min = []
    physical_max = []
    digital_min = []
    digital_max = []
    prefiltering = []
    nsr = []
    reserved2 = []
    for i in range(ns): new_ch_names.append(file.read(16).decode('ascii').strip())
    for i in range(ns): transducer.append(file.read(80).decode('ascii'))
    for i in range(ns): physical_dim.append(file.read(8).decode('ascii'))
    for i in range(ns): physical_min.append(file.read(8).decode('ascii'))
    for i in range(ns): physical_max.append(file.read(8).decode('ascii')) 
    for i in range(ns): digital_min.append((file.read(8).decode('ascii')))
    for i in range(ns): digital_max.append((file.read(8).decode('ascii')))
    for i in range(ns): prefiltering.append(file.read(80).decode('ascii'))
    for i in range(ns): nsr.append(int(file.read(8).decode('ascii')))
    for i in range(ns): reserved2.append(file.read(32).decode('ascii'))   
    record["nsr"] = nsr
    
    srate = nsr[0]/rec_in_secs
    print("srate = ", srate)
    print("Channel names: ", new_ch_names)
    print("Number of channels = ", len(new_ch_names))
    print("End of header, position = ", file.tell())
   
    # Data streams
    nsr_last = nsr[ns-1] # number of annotation entries in each record
    data = []
    edf_annotation = []
    for s in range(ns):
        data.append([])
    print("n records = ", nr)
    m = sum(nsr[0:ns-1])  # add up the number of seconds per record (nsr) for each signal
    btyes_per_record = m*2  # 2 bytes per data entry; m is the number of bytes per record
    print("Reading info for %d records ... " %(nr))
    for r in range(nr):  # Loop over records
        file.seek(btyes_per_record, 1)
        if r%1000 == 0: print("... record ", r)
        #for s in range(ns-1):
        #    data[s].extend(np.fromfile(file, dtype=np.int16,sep="", count=nsr[s]))
        
        edf_annotation.append(np.fromfile(file, dtype=np.byte,sep="", count=2*nsr_last))
       
    # Let's scale the data
    scale = 1.0
    
    physical_min = np.asarray(physical_min, dtype=float)
    physical_max = np.asarray(physical_max, dtype=float)
    digital_min = np.asarray(digital_min, dtype=float)
    digital_max = np.asarray(digital_max, dtype=float)
    print("phys min/max = ", physical_min[0], physical_max[0])
    print("digi min/max = ", digital_min[0], digital_max[0])
    m = np.zeros(ns-1)
    for s in range(ns-1):
        if   'mV' in physical_dim[0]: scale = 1.0e-3
        elif 'uV' in physical_dim[0]: scale = 1.0e-6
        elif 'nV' in physical_dim[0]: scale = 1.0e-9
        m[s] = (physical_max[s]-physical_min[s])/(digital_max[s]-digital_min[s])
        b = 0.0 #physical_min[s]
        data[s] = scale*(m[s]*np.asarray(data[s], dtype=float) + b)
    print("data:     ", data[2][0:5])
    
    # Search for discontinuities
    startdate0, starttime0 = startdate, starttime
    r0 = 0
    rn = nr-1    
    if "EDF+D" in edf_type:
        a = edf_annotation    
        rec_start = []
        for i in range(nr):
            c = a[i].tobytes().decode('ascii').split('+')
            a1 = c[1].replace("\x14", "")
            a2 = a1.replace("\x00", "")
            a3 = round(float(a2),3)
            rec_start.append(a3)
            
        t0 = rec_start[r0]
        tn = rec_start[rn]
        segments = [[t0,tn]]
        isegments = [[r0,rn]]

        s = 0
        for r in range(1,nr):
            t1 = rec_start[r]
            r1 = r
            
            dt = round(t1-t0,3)
            if dt > rec_in_secs:
                print("Gap in records %d to %d of %8.3f " % (r0, r1, t1-t0))
                segments[s][1] = t0
                isegments[s][1] = r0
                segments.append([t1,tn])
                isegments.append([r1,rn])
                
                startdate, starttime = add_seconds(startdate0, starttime0, t1)
                record["start"].append([startdate, starttime])
                s += 1
            
            t0 = t1
            r0 = r1
    else:
        t0 = 0.0
        t1 = nr*nsr[0]/srate
        segments = [[t0,t1]]
        isegments = [[r0,rn]]
                        
    nsegs = len(segments)
    record["nsegments"] = nsegs
    print("Total number of segments = ", len(segments))
    for s in range(len(segments)):
        r0 = isegments[s][0]
        rn = isegments[s][1]
        record["nrecords"].append(rn-r0+1)
        print("segment ", s+1, ": ", segments[s], "; records: ", isegments[s]," (hours: ",segments[s][1]/3600.0,")")
        print("   nrecords = ", record["nrecords"][s])
    file.close() 
           
    return record


#-----------------------------------------
# Change startdate, starttime, nr in EDF header
#-----------------------------------------
def set_header(header, record):
    # startdate: 8 bytes, starting at byte 168
    # starttime: 8 bytes, starting at byte 176
    # edf_type:  44 bytes starting at byte 192
    # nr:        8 bytes, starting at byte 236
    
    # Get the needed information
    startdate = record["start"][i][0]
    starttime = record["start"][i][1]
    nr = record["nrecords"][i]
    
    # We're writing a continuous file
    edf_type = "EDF+C"
    
    # convert the integer nr to an 8-byte string
    nr8 = str("%-8d" %(nr))  
    
    # Convert the header to a character list (strings are immutable)
    header_list = list(header.decode('ascii'))  
    
    # Replace precise characters (byte locations) in the character list
    header_list[168:176] = startdate
    header_list[176:184] = starttime
    header_list[192:197] = edf_type
    header_list[236:244] = nr8
    
    # Convert back to bytes
    bheader = bytes("".join(header_list), 'utf-8')
    
    return bheader
    

#-----------------------------------------
# Write EDF file
#-----------------------------------------
def write_edf(newfilename, bfile, header, iseg, record):
    
    ns = record["ns"]  # number of signals
    nr_seg = record["nrecords"][iseg]
    
    # Open a new binary file for writing
    file = open(newfilename, "wb")
        
    # Write the new header
    new_hdr = set_header(header, record)
    file.write(new_hdr)
    
    # Write the data records; each entry is 2 byte integer
    for r in range(nr_seg):
        for s in range(ns):
            nsr = record["nsr"][s] 
            nbytes = nsr * 2
            b_data = bfile.read(nbytes)
            file.write(b_data) 
                       
    file.close()
    

#-----------------------------------------
# main driver 
#-----------------------------------------
if __name__ == "__main__":

    # Get command line input
    infilename, rw = parse_commandline(sys.argv)
    
    # Read the file and figure out how many segments
    record = read_edfD(infilename)
    nhead = record["nhead"]
    nsegs = record["nsegments"]
    print("Number of continuous segments = ", nsegs)
        
    # Open the file again, read up until the data records    
    print ("nsegs, rw = ", nsegs, rw)
    if nsegs > 1 and rw=="write":
        # Open the original EDF+D file and read the header
        bfile = open(infilename, "rb") 
        header = bfile.read(nhead)
        print("\n----------------")
        print("Writing new continuous files ...\n")
        for i in range(nsegs):
            newfilename = infilename[:-4] + "_" + str(i+1).zfill(2) +".edf"
            print(newfilename, ": startdate, starttime = ", record["start"][i])
            write_edf(newfilename, bfile, header, i, record)
        bfile.close()
    