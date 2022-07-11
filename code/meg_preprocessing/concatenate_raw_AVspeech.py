"""
====================================================================
Concatenate M/EEG data files
====================================================================

Original script by Lucy MacGregor, edited by Heidi Solberg Okland.

Concatenate the 5 blocks of raw data per subject.
Do after maxfiltering.


"""


print __doc__

# Russell's addition for running on cluster
import sys
sys.path.insert(1, '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/sklearn/')
sys.path.insert(1, '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/pysurfer/')
sys.path.insert(1, '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/nibabel/')
sys.path.insert(1, '/imaging/local/software/mne_python/v0.14/')
sys.path.insert(1, '/imaging/local/software/freesurfer/6.0.0/')


import mne
import os.path as op
from matplotlib import pyplot as plt


# define the data path
data_path = ''

# define subject-specific subfolders
subs = []

if len(sys.argv)>1: # if in parallel mode
    print "Running subject(s) {0} now in PARALLEL mode".format(sys.argv)
    ss_idx = map(int, sys.argv[1:])
    subs_new = []
    for ii,ss in enumerate(ss_idx): # a bit cumbersome because lists cannot be used as indices
        subs_new.append(subs[ss])
    subs = subs_new
else:
    print "Running now in SERIAL mode"

# loop over subjects to read in all the blocks and concatenate them
for ss in subs:

    print "###\nReading in data files for {}...\n###".format(ss[0])

    fname_raw = data_path + ss[0] + '/' + ss[1] + '/' + 'block1_raw_ds.fif'
    fname_raw2 = data_path + ss[0] + '/' + ss[1] + '/' + 'block2_raw_ds.fif'
    fname_raw3 = data_path + ss[0] + '/' + ss[1] + '/' + 'block3_raw_ds.fif'
    fname_raw4 = data_path + ss[0] + '/' + ss[1] + '/' + 'block4_raw_ds.fif'
    fname_raw5 = data_path + ss[0] + '/' + ss[1] + '/' + 'block5_raw_ds.fif'

    raw = mne.io.read_raw_fif (fname_raw, preload=True)
    raw2 = mne.io.read_raw_fif (fname_raw2, preload=True)
    raw3 = mne.io.read_raw_fif (fname_raw3, preload=True)
    raw4 = mne.io.read_raw_fif (fname_raw4, preload=True)
    raw5 = mne.io.read_raw_fif (fname_raw5, preload=True)

    print "###\nConcatenating and saving {}...\n###".format(ss[0])

    raw = mne.concatenate_raws([raw, raw2, raw3, raw4, raw5], preload=True)

    # save (NB: overwrites if the file already exists!)
    raw.save(data_path + ss[0] + '/' + ss[1] + '/' + '/concatenated_raw.fif', overwrite=True)
