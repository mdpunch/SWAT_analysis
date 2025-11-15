# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python (ctapipe_0.24)
#     language: python
#     name: ctapipe_0.24
# ---

# %% [markdown]
# # Look into SWAT test files

# %%
# 51 is all a's
# 61,71,81 is all a's, but high rate and few?

import sys
if  not hasattr(sys, 'ps1'):

    import argparse
    parser = argparse.ArgumentParser("Read the SWAT file and produce outputs")
    
    parser.add_argument('-o', '--obs-id', help="Observation ID to read")
    parser.add_argument('-c', '--calib', default=0, type=int, help="Apply calibration pickle file with name {calib}_calibs.pickle\n"
                                                                   "  Default 0, to make calibration file from events with all telescopes hit")
    parser.add_argument('-d', '--directory', default="./triggers/2025/06/0*/*", help="Directory in which to look for files (usually of form triggers/yyyy/mm/dd/)\n"
                                                                                     "  Uses glob, so *s accepted, if argument in quotes")
    
    args = parser.parse_args()

    if args.obs_id:
        obsid =args.obs_id
    else:
        print("Need obsid!!!")
        sys.exit()
    
    calib = args.calib
    directory = args.directory

else:
    obsid = 92 # 1301 # 97 # 92 # 1003 # 1004 # 1005
    calib = 61
    directory = "./triggers/2025/06/0*/*"

# %%
import pickle
import pprint
if calib>0:
    with open(f"outputs/{calib}_calibs.pickle", 'rb') as handle:
        cal = pickle.load(handle)
    pprint.pprint(cal)

# %%
import numpy as np

# %%
from glob import glob
import re

# %%
from protozfits import File

# %%
files = glob(directory)
files.sort()
len(files),files[:10],files[-10:]

# %% [markdown]
# ## Find the ObsIDs

# %%
re.findall(".+OBSID([0-9]+).+",files[0])

pattern = re.compile(r".+OBSID([0-9]+).+")

obsids = [int(pattern.findall(file)[0]) for file in files]
obsids = sorted(list(set(obsids)))
print(f"ObsIDs available in directory {directory}:\n",obsids)


# %% [markdown]
# ## Pick an ObsID

# %%
def get_files_for_obsid(obsid):
    #for file in files:
    #    print(file+"\n" if (re.findall(f".+OBSID(?:0+){obsid}_.+",file) != []) else "",end="")

    obsid_files = [file for file in files if (re.findall(f".+OBSID(?:0+){obsid}_.+",file) != []) ]
    #obsid_files
    
    ## Separate out the `SUBARRAY` files
    
    obsid_subarray = [file for file in obsid_files if "SUBARRAY" in file]
    print(f"{len(obsid_subarray):5d} subarray files")
    
    obsid_tel = sorted(list(set(obsid_files)-set(obsid_subarray)))
    print(f"{len(obsid_tel):} telescope files")

    return obsid_subarray,obsid_tel


# %% [markdown]
# ## Look in the files

# %%
obsid_subarray,obsid_tel = get_files_for_obsid(obsid)

# %%
f_tel = File(obsid_tel[0])
f_sub = File(obsid_subarray[0])

# %%
len(f_sub.SubarrayEvents),len(f_tel.Triggers)

# %%

# %%
[(evt.event_id,evt.event_time_s,evt.event_time_qns//4,
  np.array((evt.tel_ids,evt.trigger_ids)).T) for evt in list(f_sub.SubarrayEvents[5931:5951])]

# %%
[(20000+i,tel.tel_id,tel.trigger_id,tel.trigger_time_s,tel.trigger_time_qns//4) for i,tel in enumerate(list(f_tel.Triggers[20000:20100]))]

# %%
f_sub.SubarrayEvents[1],list(f_tel.Triggers[0:1010])

# %%
if  hasattr(sys, 'ps1'):
    # Looking at outputs, SWAT takes some time to connect to all...
    f_tel = File(obsid_tel[0])
    f_sub = File(obsid_subarray[0])
    prev_sub_tel_ids = []
    prev_s,prev_qns = 0,0
    for sub in f_sub.SubarrayEvents[:10000]:
        if (len(sub.tel_ids) != len(prev_sub_tel_ids) or
            set(sub.tel_ids) != set(prev_sub_tel_ids)):
            print(sub.event_id,"Subarray tel_ids",sub.event_time_s,sub.event_time_qns,"time",sub.tel_ids)
            differ = set(sub.tel_ids.astype(int))-set(prev_sub_tel_ids)
            difftimes = (sub.event_time_s-prev_s)*1_000_000_000 + (sub.event_time_qns-prev_qns)//4
            print("  event difference", np.array(list(differ)),"times",difftimes)
            prev_sub_tel_ids = sub.tel_ids.astype(int)
        difftimes = (sub.event_time_s-prev_s)*1_000_000_000 + (sub.event_time_qns-prev_qns)//4
        print("  time difference",difftimes)
        prev_s,prev_qns = sub.event_time_s,sub.event_time_qns


# %%
from sys import maxsize
maxsize/1e9/(24*60*60)/365.25
maxsize

# %% [markdown]
# ## Loop through the files and produce outputs

# %%
ord_a = ord("a")

# %%
obsids

# %%
# Order
# Channels vs key
#   0,  1,  2,  3,  4,  5,  6,  7
#  27, 25, 23, 22, 20, 21, 26, 24

key_to_order = { 27:0, 25:1, 23:2, 22:3, 20:4, 21:5, 26:6, 24:7 }


# MP revised version with test remotely, 20250701 ... not correct for files taken!!!
# Maybe the IPs are assigned in random order after boot ?
#key_to_order = { 22:0, 27:1, 21:2, 24:3, 20:4, 25:5, 26:6, 23:7 }


# %%
import gzip
file_times = gzip.open(f"outputs/{obsid}_times.txt.gz", "wt")
file_strings = gzip.open(f"outputs/{obsid}_strings_calib_{calib}.txt.gz", "wt")

obsid_subarray,obsid_tel = get_files_for_obsid(obsid)

f_tel = File(ft:=obsid_tel.pop(0))
f_sub = File(fs:=obsid_subarray.pop(0))
print(ft,fs)

# 2.5ns spacing between deltaTs of TATS
t_step = 2.5

# In the SubarrayEvents, there is the trigger_ids, tel_ids.
#  Run through those for an event, then go through the tels, to get the times.

n_line = 0
sub_gen = f_sub.SubarrayEvents
tel_gen = f_tel.Triggers

# Start with first telescope record
tel = next(tel_gen)
t0 = tel.trigger_time_s-1_000_000_000 # Make it a second earlier for margin
sub_prev_trigger_time = t0

tel_files = {}
t_prevs = {}
e_prevs = {}

n_calibs = 0
calibs = {}

#for sub in f_sub.SubarrayEvents:
while True:
    try:
        sub = next(sub_gen)
    except StopIteration:
        f_sub.close()
        if len(obsid_subarray) == 0:
            break # No more subarray files
        next_subarray = obsid_subarray.pop(0)
        print("Nextsubarray file:",next_subarray)
        f_sub = File(next_subarray)
        sub_gen = f_sub.SubarrayEvents
        sub = next(sub_gen)        
    
    #print(sub)
    sub_tel_ids=sub.tel_ids
    sub_trigger_ids = sub.trigger_ids
    sub_trigger_time = (sub.event_time_s-t0)*1_000_000_000 + sub.event_time_qns//4
    #print("Subarray tel_ids",sub_tel_ids)
    #print("Subarray trigger_ids",sub_trigger_ids)
    sub_times={}

    tels = np.stack((sub_tel_ids,sub_trigger_ids)).T.tolist()

    #print("tels",tels)
    
    while len(tels):
        tel_id = tel.tel_id
        tel_idx = [tel.tel_id,tel.trigger_id]

    
        #print("***",tel_idx)
    
        if tel_idx in tels:
            #print("found tel_idx",tel_idx)
            
            # Open file if not already open
            if tel_id not in tel_files:
                tel_files[tel_id] = gzip.open(f"outputs/{obsid}_tel_{tel_id}.intervals.gz","wt")
                t_prevs[tel_id] = 0
                e_prevs[tel_id] = 0

            t_now_ns = tel.trigger_time_s*1_000_000_000 + tel.trigger_time_qns//4
            e_now = tel.trigger_id
            busy_status = not tel.data_available

            delta_ns = t_now_ns - t_prevs[tel_id]
            delta_evt = e_now - e_prevs[tel_id]

            tel_files[tel_id].write(f"{delta_ns:10} {delta_evt:10} {"1" if busy_status else "0"} {e_now:20} {t_now_ns:20}\n")

            t_prevs[tel_id] = t_now_ns
            e_prevs[tel_id] = e_now
            
            tels.remove(tel_idx)
            sub_times[tel.tel_id] = (tel.trigger_time_s-t0)*4_000_000_000+tel.trigger_time_qns

            continue

        # Could have a tel which has an event number higher than the subarray
        # => need to skip to the next subarray???
        # But without moving forward in the telescope file
        if tel.tel_id in sub.tel_ids:
            if tel.trigger_id > sub_trigger_ids[np.argwhere(sub_tel_ids==tel.tel_id)[0,0]]:
                #print(20*"v"+" TRIGGER AFTER SUB-ARRAY ??? "+20*"v")
                #file_strings.write(20*"v"+" TRIGGER AFTER SUB-ARRAY ??? "+20*"v"+"\n")
                continue
        else:
            #print(20*">"+" NON-COINCIDENT EVENT ??? "+20*">")
            #file_strings.write(20*">"+" NON-COINCIDENT EVENT ??? "+20*">"+"\n")
            # Are these at the same time as a missing packet?
            pass

        #print(">>>",tel_idx,tel.trigger_time_s,tel.trigger_time_qns)

        try:
            tel = next(tel_gen)
        except StopIteration:
            f_tel.close()
            if len(obsid_tel) == 0:
                break # No more telescope files
            next_tel = obsid_tel.pop(0)
            print("Next telescope file:",next_tel)
            f_tel = File(next_tel)
            tel_gen = f_tel.Triggers
            tel = next(tel_gen)

    # Sort by channel
    sub_times = dict(sorted(sub_times.items()))
    
    
    # get differences and convert from quarter_ns to ns
    min_times = min(sub_times.values())
    for key in sub_times.keys():
        sub_times[key] -= min_times
        sub_times[key] //= 4
    #print(sub_times)

    min_chan = min(sub_times, key=sub_times.get)

    if calib>0:
        #print("<<<")
        #pprint.pprint(sub_times)
        # Run through excluding 1st channel (but that shouldn't make a difference if included)
        for key in list(sub_times.keys())[1:]:
            sub_times[key] -= cal[key] - cal[min_chan] 
        #print(">>>")
        #pprint.pprint(sub_times)
 
    # re-do get differences
    min_times = min(sub_times.values())
    for key in sub_times.keys():
        sub_times[key] -= min_times

    #if len(sub_times)==8:
    #    break

    if calib>0:
        #print("===")
        #pprint.pprint(sub_times)
        pass

    
    '''
    if max(sub_times.values())>9:
        print("... FAR DISTANT DELTA_T ...")
        print(sub_times)
    '''

    times_str = "........"
    lts = list(times_str)
    for key in sub_times.keys():
        lts[key_to_order[key]] = f"{sub_times[key]}"[-1]    
    times_str = "".join(lts)
    
    #print(trig_str)
    file_times.write(times_str+"\n")

    '''
    #chr(ord_a + np.array(list(sub_times.values()))//4)
    # +np.round((t4th_ns-min_time)/cnt_per_ns/clk_period).astype('int')
    trig_str = "........"
    lts = list(trig_str)
    for key in sub_times.keys():
        lts[key_to_order[key]] = chr(ord_a + np.floor_divide(sub_times[key],t_step).astype('int'))  
    trig_str = "".join(lts)

    # Now, check if there is a better solution, by looking for minimum
    nearest_grid_times = {}
    sub_times_values = np.array(list(sub_times.values()))
    rmsq_min = np.inf
    for delta in np.linspace(-t_step/2,t_step/2,21):
        grid = np.round((sub_times_values+delta)/t_step).astype('int')
        grid_diff = sub_times_values-grid*t_step
        rmsq = np.sqrt(np.sum(grid_diff**2))
        if rmsq < rmsq_min:
            rmsq_min, delta_min, grid_min = rmsq, delta, grid
        elif rmsq == rmsq_min:
            delta_min2 = delta
    trig_str_min = "........"
    ltsm = list(trig_str_min)
    stv_keys = list(sub_times.keys())
    for key,gm in zip(stv_keys,grid_min):
        ltsm[key_to_order[key]] = chr(ord('a')+gm)
    trig_str_min = "".join(ltsm)

    if trig_str != trig_str_min:
        if n_line<10000 or n_line%100_000==0:
            if delta_min > 0 or delta_min2 < 0:
                print("        >> minimum better than nearest:",trig_str_min,trig_str,"minimum from",delta_min,"--",delta_min2)
        trig_str = trig_str_min

    '''
    
    # Rounding tested to be the best solution
    trig_str = "........"
    lts = list(trig_str)
    for key in sub_times.keys():
        lts[key_to_order[key]] = chr(ord('a')+
                                     np.round((sub_times[key])/t_step).astype('int'))
    trig_str = "".join(lts)


    dt = sub_trigger_time - sub_prev_trigger_time
    #print(trig_str)
    file_strings.write(trig_str+f" {dt}\n")
    sub_prev_trigger_time = sub_trigger_time

    # Try making the calibration w.r.t. 0th channel
    zeroth_chan = list(sub_times.keys())[0]
    if calib==0:
        if len(sub_times)==8:
            for key in sub_times.keys():
                calibs[key] = calibs.get(key, 0) + sub_times[key]-sub_times[zeroth_chan]
            n_calibs += 1
    
    if n_line<10000:
        print(trig_str)
    elif n_line%100_000==0:
        print(trig_str,n_line)
    elif n_line>1_000_000_000:
        break
    #print(80*"-")
    n_line += 1

# %%
obsid

# %%
sub_times[key]

# %%
key_to_order[key]

# %%
sub_times_values_reordered

# %%
if calib==0:
    cal = {k:v/n_calibs for k,v in calibs.items()}
    print(n_calibs,cal)
    with open(f"outputs/{obsid}_calibs.pickle", 'wb') as handle:
        pickle.dump(cal, handle, protocol=pickle.HIGHEST_PROTOCOL)

# %%
sub_times

# %%
key_to_order

# %%
