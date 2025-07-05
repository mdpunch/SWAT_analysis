#! /usr/bin/env python3.7
# This reads in the data file produced by just one card, then makes interval distribution
# e.g. "./Make_Intervals.py 55000stamps.bin 100 | gzip -v >  Intervals.txt.gz "

import gzip
import sys

import struct


def unpackdata(data,chan):
    '''
    Input:
            data, the data block corresponding to an event
    Output:
            time, the time of the current event
    '''

    ''' The data should be this structure of data 
    {
        tUInt64 uctsTimeStamp;
        tUInt32 uctsAddress;               /// Note, in the unpack, this takes 4 values xxx.yyy.zzz.ddd
        tUInt32 eventCounter;
        tUInt32 busyCounter;
        tUInt32 ppsCounter;
        tUInt32 clockCounter;
        tUInt8 triggerType;                 // For TIB, this is the first 7 bits of SPI
        tUInt8 whiteRabbitResetBusyStatus;  // A few status flags here, see below
        tUInt8 stereoPattern;               // For TIB, this is the next 8 bits of the SPI string
        tUInt8 numInBunch;                  // This information only needed for debugging
        tUInt32 cdtsVersion;                // Should be x.y.z where x is 2bytes and y, z are a byte each
    } __attribute__((__packed__)) tExtDevicesCDTSMessage;
    '''
    try:
      # 4 + 4x4 + 1x8 + 4x1 = 
      tuple_of_data = struct.unpack("QBBBBIIIIBBBBI", data)
    except struct.error:
      #print("Struct Error")
      return -99,-1,-1,-1#,tuple_of_data
    #print("*",tuple_of_data)
    #print("**",tuple_of_data[4])
    read_chan = tuple_of_data[4]

    #print(tuple_of_data)

    time = tuple_of_data[0]
    ev_num = tuple_of_data[5]
    bz_num = tuple_of_data[6]
    pps_ct = tuple_of_data[7]
    clk_ct = tuple_of_data[8]
    trg_typ = tuple_of_data[9]
    str_pat = tuple_of_data[11]
    num_bch = tuple_of_data[12]

    if chan == 0 or chan == read_chan:
      return time,ev_num,bz_num,pps_ct,clk_ct,trg_typ,str_pat,num_bch#,tuple_of_data
    return -1,-1,-1,-1,-1,-1,-1,-1#,tuple_of_data


if len(sys.argv) < 2:
    print("ERROR: You must give the binary file with the timestamps as an argument!")
    exit()

fname = sys.argv[1]    

chan = 0 
if len(sys.argv) <3:
    print("Assuming time stamps in file are all from the same TiCkS board")
else:
    chan = int(sys.argv[2])

# Modified to assume data in ../TiCkS_data by default
try:
    #f = open(fname)
    #txt = f.readlines()
    fstamps = open("../../TiCkS/TiCkS_data/"+fname,"rb")
    #fstamps = open(fname,"rb")
except IOError:
    try:
        fstamps = open(fname,"rb")
    except IOError:
        print("ERROR: File ","../../TiCkS/TiCkS_data/",fname,"not found, nor ",fname,"!")
        #print("ERROR: File ",fname,"not found!")
        exit()

event_size = 36
data = fstamps.read(event_size)
tprev, eprev, bprev, pprev, cprev, trg_typ_prev, str_pat_prev, num_bch = unpackdata(data,0)
# tats event #, MSB is reserved for ????
tats_prev = ((str_pat_prev&127)<<7) + (trg_typ_prev)
#tats_prev = (str_pat_prev<<7) + (trg_typ_prev)
tnow = 0 
while(tnow!=-99):

    data = fstamps.read(event_size)
    tnow, enow, bnow, ppsnow, clknow, trg_typ, str_pat, num_bch = unpackdata(data,chan)
    #tnow, enow, bnow, tuple_of_data = unpackdata(data,chan)
    if tnow == -1: # Time not from the channel we want
        continue
    Interval = tnow-tprev
    #print(Interval,enow-eprev,bnow-bprev,num_bch)
    #print(f'{Interval:10} {enow-eprev:1} {bnow-bprev:1}')
    tats = ((str_pat&127)<<7) + (trg_typ)
    #tats = (str_pat<<7) + (trg_typ)
    print(f'{Interval:10} {enow+bnow-eprev-bprev:1} {enow-eprev:1} {bnow-bprev:1} '
          f'{(tats-tats_prev)%16384} {tnow:40} {enow:10} {bnow:10} {tats}',end='') #16384 is 2^14, 32768 is 2^15
    #print(">>> clk vs. clk part of ts, diff:",clknow,"vs.",int((tnow%1000000000)/100),", clk-clk_ts =", 
    #          clknow-int((tnow%1000000000)/100)) 
    clk_diff = clknow-int((tnow%1000000000)/100)
    if clk_diff >= 9999998 and ppsnow != pprev:  # If roll-over of PPS, but clock not rolled-over yet
      clk_diff -= 10000000
    print(f"   {ppsnow:20} {clknow:10} {clk_diff:10}") 

    #print("**",tats,str_pat,trg_typ)

    #if enow-eprev>1:
    #    print(tuple_of_data)
    tprev,eprev,bprev = tnow,enow,bnow
    trg_typ_prev,str_pat_prev = trg_typ,str_pat
    pprev,cprev = ppsnow,clknow
    tats_prev = tats


fstamps.close()
