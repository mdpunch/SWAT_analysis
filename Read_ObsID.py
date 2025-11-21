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
import ast

default_keyorder = { 27:0, 25:1, 23:2, 22:3, 20:4, 21:5, 26:6, 24:7 }

if  not hasattr(sys, 'ps1'):

    import argparse
    parser = argparse.ArgumentParser("Read the SWAT file and produce outputs")
    
    parser.add_argument('-o', '--obs-id', help="Observation ID to read")
    parser.add_argument('-c', '--calib', default=0, type=int, help="Apply calibration pickle file with name {calib}_calibs.pickle.\n"
                                                                   "  Default 0, to make calibration file from events with all telescopes hit.")
    parser.add_argument('-d', '--directory', default="./triggers/2025/06/0*/*", help="Directory in which to look for files (usually of form triggers/yyyy/mm/dd/).\n"
                                                                                     "  Uses glob, so *s accepted, if argument in quotes.")
    parser.add_argument('-k', '--keyorder', default="",help=f"Dictionary (as string) of the order of the tels in the TATS string.  e.g. {str(default_keyorder)}. Looks in file \"key_to_order.txt\" if not, or defaults this.")
    
    args = parser.parse_args()

    if args.obs_id:
        obsid =args.obs_id
    else:
        print("Need obsid!!!")
        sys.exit()
    
    calib = args.calib
    directory = args.directory
    if len(args.keyorder) == 0: # No keyprder provided on command line, try for file, and if no file default to default
        try:
            with open("key_to_order.txt", "r") as key_to_order_file:
                file_key_to_order = key_to_order_file.readline()
        except:
            print(f"No keyOrder file \"key_to_order.txt\", using default {str(default_keyorder)}.")
            key_to_order = default_keyorder
        else:
            try:
                key_to_order = ast.literal_eval(file_key_to_order)
            except:
                print("Error interpreting keyOrder string from \"key_to_order.txt\" as a dictionary:",file_key_to_order)
                print("Exiting...")
                sys.exit()
            print(f"Using keyOrder string from \"key_to_order.txt\" as a dictionary: {key_to_order}")

    else:
        try:
            key_to_order = ast.literal_eval(args.keyorder)
        except:
            print("Error interpreting keyOrder string argument as a dictionary:",args.keyorder)
            print("Exiting...")
            sys.exit()
    

else:
    obsid = 92 # 1301 # 97 # 92 # 1003 # 1004 # 1005
    calib = 61
    directory = "./triggers/2025/06/0*/*"
    key_to_order = default_keyorder

if len(key_to_order) != 8:
    print("Not enough channels in keyOrder dictionary",key_to_order)
    print("Exiting...")
    sys.exit()

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
print(f"Number of files found in directory {directory}:",len(files))
files[:10],files[-10:]

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
  np.array((evt.tel_ids_with_trigger,evt.trigger_ids)).T) for evt in list(f_sub.SubarrayEvents[5931:5951])]

# %%
[(20000+i,tel.tel_id,tel.trigger_id,tel.trigger_time_s,tel.trigger_time_qns//4) for i,tel in enumerate(list(f_tel.Triggers[20000:20100]))]

# %%
f_sub.SubarrayEvents[1],list(f_tel.Triggers[0:1010])

# %%
f_sub.SubarrayEvents[1].trigger_ids

# %%
if  hasattr(sys, 'ps1'):
    # Looking at outputs, SWAT takes some time to connect to all...
    f_tel = File(obsid_tel[0])
    f_sub = File(obsid_subarray[0])
    prev_sub_tel_ids_with_trigger = []
    prev_s,prev_qns = 0,0
    for sub in f_sub.SubarrayEvents[:10000]:
        if (len(sub.tel_ids_with_trigger) != len(prev_sub_tel_ids_with_trigger) or
            set(sub.tel_ids_with_trigger) != set(prev_sub_tel_ids_with_trigger)):
            print(sub.event_id,"Subarray tel_ids",sub.event_time_s,sub.event_time_qns,"time",sub.tel_ids_with_trigger)
            differ = set(sub.tel_ids_with_trigger.astype(int))-set(prev_sub_tel_ids_with_trigger)
            difftimes = (sub.event_time_s-prev_s)*1_000_000_000 + (sub.event_time_qns-prev_qns)//4
            print("  event difference", np.array(list(differ)),"times",difftimes)
            prev_sub_tel_ids_with_trigger = sub.tel_ids_with_trigger.astype(int)
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

#key_to_order = { 27:0, 25:1, 23:2, 22:3, 20:4, 21:5, 26:6, 24:7 }


# MP revised version with test remotely, 20250701 ... not correct for files taken!!!
# Maybe the IPs are assigned in random order after boot ?
#key_to_order = { 22:0, 27:1, 21:2, 24:3, 20:4, 25:5, 26:6, 23:7 }


# %%
## DEPRECATED ... DRY
def get_next_tel_trigger(obsid_tel):
    ''' get next telescope from list of triggers files

    obsid_tel: list of telescope trigger files
        
    returns tel
    '''

    
    # Loop over files
    while True:
        
        try:
            ft = obsid_tel.pop(0)
        except IndexError:
            # No more telescope files
            return # bare return equivalent to StopIteration
        
        f_tel = File(ft)

        tel_gen = f_tel.Triggers
        
        # Loop over triggers in files
        while True:
            try:
                tel = next(tel_gen)
                yield tel
            except StopIteration:
                f_tel.close()
                break 


# %%
## DEPRECATED ... DRY
def get_next_subarray_event(obsid_sub):
    ''' get next subarray events from list of subarray files

    obsid_sub: list of subarray events files
        
    returns sub
    '''

    
    # Loop over files
    while True:
        
        try:
            ft = obsid_sub.pop(0)
        except IndexError:
            # No more subarray files
            return # bare return equivalent to StopIteration
        
        f_sub = File(fs)

        sub_gen = f_sub.SubarrayEvents
        
        # Loop over subarray events in files
        while True:
            try:
                sub = next(sub_gen)
                yield sub
            except StopIteration:
                f_sub.close()
                break


# %%
def get_next_record(filelist,thing_name):
    ''' get next thing in record from list of files

    filelist:    
        list of files
    thing_name:  
        attribute to get from the object (string)
        
    returns the thing
    '''

    
    # Loop over files
    while True:
        
        try:
            filename = filelist.pop(0)
        except IndexError:
            # No more subarray files
            return # bare return equivalent to StopIteration
        
        file_handle = File(filename)
        print("Opened",filename)

        thing_gen = getattr(file_handle,thing_name)
        
        # Loop over records in files
        while True:
            try:
                thing = next(thing_gen)
                yield thing
            except StopIteration:                
                file_handle.close()
                break


# %%
import copy

# %%
import gzip
file_times = gzip.open(f"outputs/{obsid}_times.txt.gz", "wt")
file_strings = gzip.open(f"outputs/{obsid}_strings_calib_{calib}.txt.gz", "wt")

obsid_subarray,obsid_tel = get_files_for_obsid(obsid)

#get_next_subarray_event_gen = get_next_subarray_event(obsid_subarray)
get_next_subarray_event_gen = get_next_record(obsid_subarray,"SubarrayEvents")

#get_next_tel_trigger_gen = get_next_tel_trigger(obsid_tel)
get_next_tel_trigger_gen = get_next_record(obsid_tel,"Triggers")

# 2.5ns spacing between deltaTs of TATS
t_step = 2.5

# In the SubarrayEvents, there is the trigger_ids, tel_ids.
#  Run through those for an event, then go through the tels_in_sub, to get the times.

n_line = 0

# Start with first telescope record
#tel = next(get_next_tel_trigger_gen)
t0 = 0
sub_prev_trigger_time = t0

tel_files = {}
t_prevs = {}
e_prevs = {}

n_calibs = 0
calibs = {}

#for sub in f_sub.SubarrayEvents:

# Loop over subarray events
count = 0
while True:
    
    try:
        sub = next(get_next_subarray_event_gen)
    except StopIteration:
        break

    sub_tel_ids_with_trigger=sub.tel_ids_with_trigger
    sub_trigger_ids = sub.trigger_ids
    sub_trigger_time = (sub.event_time_s-t0)*1_000_000_000 + sub.event_time_qns//4
    #print("Subarray tel_ids",sub_tel_ids_with_trigger)
    #print("Subarray trigger_ids",sub_trigger_ids)
    sub_times={}

    # Make a list of the telescopes with trigger in the subarray for this event, 
    #    and their trigger IDs for this event.
    tels_in_sub = np.stack((sub_tel_ids_with_trigger,sub_trigger_ids)).T.tolist()

    #print(len(tels_in_sub),"tels_in_sub",tels_in_sub)
    
    # Run through the telescope trigger file and 
    # search for the corresponding tel/trigger_id in the telescope triggers file
    # Fill the trigger corresponding times, until all tels in subarray found
    
    while tel := next(get_next_tel_trigger_gen):
        
        count += 1
        #if count>5000000:
        #    break
        #print(count,end=" ")
        
        # Start with the current tel/trigger_id in the telescope triggers file
        tel_id = tel.tel_id
        tel_idx = [tel.tel_id,tel.trigger_id]
    
        #print("***",tel_idx)
    
        if tel_idx in tels_in_sub:
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
            
            tels_in_sub.remove(tel_idx)
            sub_times[tel.tel_id] = (tel.trigger_time_s-t0)*4_000_000_000+tel.trigger_time_qns

            # Note: if len(tels_in_sub) is now zero, we exit the loop,
            #   treat the subarray_times and then search for the next subarray event.
            if len(tels_in_sub) == 0:
                break

            # Add a check that there isn't something crazy going on.
            # If a telescope trigger is missing, then what?

    # convert from quarter_ns to ns
    sub_times = {k:v//4 for k,v in sub_times.items()}
    # Sort by channel
    sub_times = dict(sorted(sub_times.items()))
    # Backup for debug print
    sub_times_raw = copy.copy(sub_times)
    
    
    # get differences
    min_times = min(sub_times.values())
    for key in sub_times.keys():
        sub_times[key] -= min_times
        #sub_times[key] //= 4
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

    
    #if max(sub_times.values())>9:
    #    print("... FAR DISTANT DELTA_T ...")
    #    print(sub_times)
    

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
    file_strings.write(trig_str+f" {dt:10d}    "+str(sub_times_raw)+"\n")
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
if calib==0:
    cal = {k:v/n_calibs for k,v in calibs.items()}
    print(n_calibs,cal)
    with open(f"outputs/{obsid}_calibs.pickle", 'wb') as handle:
        pickle.dump(cal, handle, protocol=pickle.HIGHEST_PROTOCOL)

# %%
