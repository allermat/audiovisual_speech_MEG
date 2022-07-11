
""" 
========================================================
Check EEG locations
========================================================

When EEG channels > 60 as at CBU, the EEG channel location obtained from the Polhemus digitiser is not copied properly to
Neuromag acquisition software. Therefore one must apply mne_check_eeg_locations to the data. Do this as early as possible in the
processing pipeline.

http://imaging.mrc-cbu.cam.ac.uk/meg/AnalyzingData/MNE_FixingFIFF

"""

import os
from os.path import isfile

path = ''

subjfolder = []

filename = 'concatenated_raw.fif'

for ss in subjfolder:
     now_subj = ' '.join(ss) # convert the list item to a string

     data_file_in = path + now_subj + filename

     check_file = isfile(data_file_in) # check if the file exists

     if check_file:
            print "###\nThe file {} exists, good!\n###".format(data_file_in)
            check_eeg = 'mne_check_eeg_locations --fix --file ' + data_file_in
            print "###\nFixing EEG locations for {}...\n###".format(data_file_in)
            os.system(check_eeg)
     else:
            print "###\nThe file {0} does NOT exist. Please check and run again.\n###".format(data_file_in)





