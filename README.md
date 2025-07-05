# SWAT files Analysis

## Read_ObsID

Python notebook or script to read files from a given ObsID.

On the command line, can run this with 
`python Read_ObsID.py {obsID}`

Reads through the Subarray and Telescope files for the ObsID, and creates:

* Telescope files, one per channel, with name `{obsid}_tel_{tel}.intervals.gz` with lines:
    - `deltaT_ns delta_evt# busy_status(0/1) event_num T_now_ns`:20}\n")`
* Coincidence time files with name `{obsid}_times.txt.gz` with lines
    - \#chan chars, with `.` if telescope not hit, and `n` for the last char of the deltaT to the 1st telescope hit
* Coincidence string files with name `{obsid}_strings_calib_{calib}.txt.gz` with lines
    - \#chan chars, with `.` if telescope not hit, and `c` for the char `ord('a')+deltaT/t_step` to the 1st telescope hit  
    where `t_step` is the step of the TATS.
* If `calib` is `0`, then it makes calibration pickle file with name `{obsid}_calibs.pickle`,   
  which contains the average deltaT from chan_0 for the other channels for events where all telescopes are hit
    - This file can be shown easily with `./show_calibs.py {obsid}_calibs.pickle`
* If `calib` is `>0`, then the corresponding `{calib}_calibs.pickle` file is applied to the times

The latter strings/times are reordered from the channels as seen by the SWAT, to be in TATS order.

This reordering is determined at runtime, by setting only one channel in TATS, and seeing which channel is seen in SWAT (or just `tcpdump`, show which channels are active).

For example, this gives:
```
# Channels vs key
#   0,  1,  2,  3,  4,  5,  6,  7
#  27, 25, 23, 22, 20, 21, 26, 24

key_to_order = { 27:0, 25:1, 23:2, 22:3, 20:4, 21:5, 26:6, 24:7 }
```

*Note:* The calibration from one obs_id was applied to another obs_id with simultaneous telescopes, to check that the calibration was working!
 

## Plot_Deltas

Make some plots of rates, deltaT distribution with exponential fit, missing events distribution (and number).  Best to start the plot at 500ns, since the fit will be disturbed by the deadtime gap.

```
usage: Plot distribution from file, mainly for Intervals [-h] [-f FILE] [-o OBS_ID] [-c CHANNEL] [-r RANGE] [-b BINS] [-l | --log | --no-log]
                                                         [-n NUM] [-t TITLE] [--cartouche CARTOUCHE] [-y YMAX]

options:
  -h, --help            show this help message and exit
  -f, --file FILE       Gzipped file to read from, Intervals.txt.gz, ignored if obs-id/channel
  -o, --obs-id OBS_ID   Observation ID of file
  -c, --channel CHANNEL
                        Telescope Channel of file
  -r, --range RANGE     Range [min, max] for plot, default from plot
  -b, --bins BINS       Number of bins, default 10 or range if -1, max 1000
  -l, --log, --no-log
  -n, --num NUM         Number of events to plot, or range, default all,
  -t, --title TITLE     Title for plot
  --cartouche CARTOUCHE
                        Legend location
  -y, --ymax YMAX       Maximum of y-axis scale
```

e.g. `python ./Plot_Deltas.py -o 91 -c 22 -n 1000000 -r [500,150000] -b 500 -n [0,100000000]`