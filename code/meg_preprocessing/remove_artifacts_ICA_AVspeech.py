"""
====================================================
Preprocessing MEG and EEG data using filters and ICA
====================================================

Written by Lucy MacGregor/based on MNE tutorials, modified by Heidi Solberg Okland

This script will find and remove eyeblink and heartbeat artefacts in your M/EEG data using ICA.
It allows you to:
* Filter MEG and EEG data and subsequently compute separate ICA solutions for MEG and EEG
* Create an html report that includes diagnostic plots of the preprocessing steps
* Specify the maximum numbers of components would wish to reject (n_max_eog and n_max_ecg)
* Plot figures showing components, topographies, reconstruction of signal for eyeblinks and heartbeats separately
* Type in which components to reject after having visually inspected the plots
* Save and apply ICA solutions (optional)

NB: Will not work on CBU-collected EEG data unless you've checked EEG electrode names using 'mne_check_eeg_locations'!


"""

print __doc__

# Russell's addition for running on cluster
import sys
sys.path.insert(1, '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/')
sys.path.insert(1, '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/sklearn/')
sys.path.insert(1, '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/pysurfer/')
sys.path.insert(1, '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/nibabel/')
sys.path.insert(1, '/imaging/local/software/mne_python/v0.11/')

import os
import numpy as np
import matplotlib.pyplot as plt

import mne
from mne.io import Raw
from mne.preprocessing import create_eog_epochs, create_ecg_epochs
from mne.preprocessing.ica import ICA, run_ica
 
mne.utils.set_log_level('WARNING')

plt.ion() # interactive plotting

whoami = os.environ['USER']
data_path = ''
fig_path = ''

# specify ichannel selection and whether or not to do ecg and/or eog
channels = ['meg', 'eeg']
print "Doing both EEG and MEG."
#channels = ['meg']
#print "Doing MEG only."
do_ecg = 'yes' 
do_eog = 'yes'

# define subject-specific subfolders
subs = []

# -----------> NB: remember to specify bad EEG electrodes down below in the loop!

# run in parallel or serial mode?
if len(sys.argv)>1: # if in parallel mode
    print "###\nRunning subject(s) {0} now in PARALLEL mode...\n###".format(sys.argv)
    ss_idx = map(int, sys.argv[1:])
    subs_new = []
    for ii,ss in enumerate(ss_idx): # a bit cumbersome because lists cannot be used as indices
        subs_new.append(subs[ss])
    subs = subs_new
else:
    print "###\nRunning now in SERIAL mode...\n###"

n_subjects = len(subs)


# loop over subjects
for ss in subs: 

  ####################### Read in and filter raw data #######################

  print "###\nPreparing subject " + ss[0] + "...\n###"

  fname_raw_in = data_path + ss[0] + '/' + ss[1] + '/concatenated_raw.fif'
  fname_raw_out = data_path + ss[0] + '/' + ss[1] + '/concatenated_icaed_raw.fif'
  fname_ica_out = data_path + ss[0] + '/' + ss[1] + '/ica.fif'   # file for ICA decomposition
  subID = ss[0]

  print 'Raw file is {}' .format(fname_raw_in)
  print 'ICA file to be saved is {}\n\n' .format(fname_ica_out)

  ## read raw data
  raw = Raw(fname_raw_in, preload=True)

  # filter raw data
  raw.filter(1, 45, n_jobs=1)

  # specify bad channels (only EEG ones neccessary - Maxfilter should have taken care of bad MEG channels)

  if ss[0]=='meg17_0066':
    raw.info['bads'] = ['EEG003', 'EEG030', 'EEG009', 'EEG070', 'EEG073'] 

  if ss[0]=='meg17_0128':
    raw.info['bads'] = ['EEG036', 'EEG037', 'EEG039', 'EEG050', 'EEG065','EEG068']

  if ss[0]=='meg17_0161':
    raw.info['bads'] = ['EEG006', 'EEG009', 'EEG019', 'EEG026', 'EEG034', 'EEG044', 'EEG044', 'EEG068', 'EEG070']

  if ss[0]=='meg17_0166':
    raw.info['bads'] = ['EEG002', 'EEG006', 'EEG018', 'EEG032', 'EEG033', 'EEG040', 'EEG073']

  if ss[0]=='meg17_0169':
    raw.info['bads'] = ['EEG009', 'EEG018', 'EEG028', 'EEG029', 'EEG031', 'EEG034', 'EEG073']

  if ss[0]=='meg17_0170':
    raw.info['bads'] = ['EEG002', 'EEG029', 'EEG032', 'EEG033', 'EEG040', 'EEG050', 'EEG051', 'EEG052', 'EEG065', 'EEG073']

  if ss[0]=='meg17_0176':
    raw.info['bads'] = ['EEG006', 'EEG073']

  if ss[0]=='meg17_0178':
    raw.info['bads'] = ['EEG011', 'EEG038', 'EEG073']

  if ss[0]=='meg17_0182':
    raw.info['bads'] = ['EEG002', 'EEG033', 'EEG037', 'EEG043', 'EEG065', 'EEG068', 'EEG073']

  if ss[0]=='meg17_0191':
    raw.info['bads'] = ['EEG005', 'EEG010', 'EEG018', 'EEG073']

  if ss[0]=='meg17_0197':
    raw.info['bads'] = ['EEG003', 'EEG004']

  if ss[0]=='meg17_0201':
    raw.info['bads'] = ['EEG002', 'EEG032', 'EEG033', 'EEG045', 'EEG073']

  if ss[0]=='meg17_0202':
    raw.info['bads'] = ['EEG009', 'EEG029', 'EEG039', 'EEG070', 'EEG073']

  if ss[0]=='meg17_0203':
    raw.info['bads'] = ['EEG002', 'EEG032', 'EEG033', 'EEG036', 'EEG051', 'EEG073']

  if ss[0]=='meg17_0214':
    raw.info['bads'] = ['EEG007', 'EEG033']


  print "You have defined the following bad electrodes for {}:\n{}" .format(ss[0], raw.info['bads'])


  ####################### ICA decomposition #######################

  # choose channel types for the ICA decomposition (M/EEG only)
  picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False, ecg=False,
                  stim=False, exclude='bads') # eeg = True
  #picks = mne.pick_types(raw.info, meg=True, eeg=False, eog=False, ecg=False,
  #                           stim=False, exclude='bads') # eeg = False

  # prepare ICA decomposition
  ica = ICA(n_components=0.90, n_pca_components=None, max_pca_components=None,
    noise_cov=None, random_state=0)
  print ica

  # read whole raw data file
  start, stop = None, None

  # set high rejection parameters to avoid computing ICA on too artifacts segments
  reject = dict(mag=4e-12, grad=4000e-13, eeg=200e-6)
  #reject = dict(mag=4e-12, grad=4000e-13)

  # decompose sources for raw data
  ica.fit(raw, start=start, stop=stop, picks=picks, decim=None, reject=reject)
  print ica

  # estimate ICA sources given unmixing matrix
  sources = ica.get_sources(raw, start=start, stop=stop)

  # setup reasonable time window for inspection, get indices for time samples
  #start_plot, stop_plot = raw.time_as_index([1, 5])
  start_plot, stop_plot = None, None

  # plot ICA component time courses (note: 'picks' here represent ICA components, not channels)
  fig = ica.plot_sources(raw, picks=None, start=start_plot, stop=stop_plot);
  plt.savefig(fig_path + subID + '_ICA_components')

  # define max number of components to reject (you may wish to edit this!)
  n_max_eog, n_max_ecg = 2, 3

  ####################### Find and remove eyeblinks #######################

  if do_eog =='yes':
    print "###\nDoing EOG.\n###"
    eog_tmin, eog_tmax = -0.5, 0.5 # how much time do we want around the blinks to make sure we capture them?
    reject_eog = {'mag': 5e-12, 'grad': 5000e-13}

    print "###\nUsing EOG062 as EOG channel.\nEpoching around the EOG peaks with tmin = {} and tmax = {}.\n###" .format(eog_tmin, eog_tmax) 
    # find EOG peaks and create epochs around them according to eog_tmin and eog_tmax
    eog_epochs = create_eog_epochs(raw, ch_name='EOG062', tmin=eog_tmin, tmax=eog_tmax,
	                            l_freq=1, h_freq=10, reject=reject_eog)
    # find ICA components that correlate highly with EOG channel using epochs or raw data
    eog_inds, eog_scores = ica.find_bads_eog(eog_epochs, ch_name='EOG062', threshold=3.0) # 3.0 is std. reduce threshold to find more components
  
    #eog_inds, eog_scores = ica.find_bads_eog(raw, l_freq=1, h_freq=10)
    eog_inds[:n_max_eog]

    # or manually specify the ICA components that correlate with the EOG
    # eog_inds = [ 3]

    if len(eog_inds)>0: #only if EOG component has been identified
      # plot scores of ICA components, mark the ones to be excluded
      fig = ica.plot_scores(eog_scores, exclude=eog_inds, labels='eog',
                                show=True, title='EOG')
      plt.savefig(fig_path + subID + '_eog_components')


      # plot topographies of ICA components to be excluded (mag and grad separately
      fig = ica.plot_components(eog_inds, ch_type='mag',
                                      title='MAG', colorbar=True, show=True)
      plt.savefig(fig_path + subID + '_eog_topo_mag')


      fig = ica.plot_components(eog_inds, ch_type='grad',
                                      title='GRAD', colorbar=True, show=True)
      plt.savefig(fig_path + subID + '_eog_topo_grad')


      fig = ica.plot_components(eog_inds, ch_type='eeg',
                                      title='EEG', colorbar=True, show=True)
      plt.savefig(fig_path + subID + '_eog_topo_eeg')

      answer = raw_input("Do you want to exclude EOG component(s) at this time? (y/n)\n")

      if answer == 'y': # yes, we want to go ahead!
        # choose which components to mark as eyeblink artefacts
        eog_reject = raw_input('\nType the indices for the eyeblink component(s) you want to reject separated by commas (e.g. 0, 3):\n')
        eog_reject = map(int,eog_reject.split(',')) # convert string input to a list of integers
        # put choice(s) in exclusion list 
        ica.exclude.extend(eog_reject)
      else:
        print "OK, not excluding any EOG components."
        eog_reject = None

    else:
      print "No EOG components found! Maybe you should try tweaking the parameters."
      eog_reject = None

  else:
    print "###\nNOT doing EOG.\n###"

  ####################### Find and remove heartbeat #######################

  if do_ecg=='yes':
    print "###\nDoing ECG.\n###"
    ecg_tmin, ecg_tmax = -0.5, 0.5
    reject_ecg = {'mag': 5e-12, 'grad': 5000e-13}
    ecg_epochs = create_ecg_epochs(raw, ch_name='ECG063', tmin=ecg_tmin, tmax=ecg_tmax,
                                     keep_ecg=False, reject=reject_ecg)
    
    # find ICA components that correlate highly with ECG channel using epochs or raw data
    #ecg_inds, ecg_scores = ica.find_bads_ecg(ecg_epochs)
    ecg_inds, ecg_scores = ica.find_bads_ecg(ecg_epochs, ch_name = 'ECG063', method= 'ctps')
    ecg_inds = ecg_inds[:n_max_ecg]

    if len(ecg_inds)>0: # only if ECG component has been identified
      # plot scores of ICA components, mark the ones to be excluded
      fig = ica.plot_scores(ecg_scores, exclude=ecg_inds, labels='ecg',
                                show=True, title='ECG')
      plt.savefig(fig_path + subID + '_ecg_components')

      # plot topographies of ICA components to be excluded (mag and grad separately)
      fig = ica.plot_components(ecg_inds, ch_type='mag',
                                    title='MAG', colorbar=True, show=True)
      plt.savefig(fig_path + subID + '_ecg_topo_mag')

      fig = ica.plot_components(ecg_inds, ch_type='grad',
                                    title='GRAD', colorbar=True, show=True)
      plt.savefig(fig_path + subID + '_ecg_topo_grad')


      fig = ica.plot_components(ecg_inds, ch_type='eeg',
                                   title='EEG', colorbar=True, show=True)
      plt.savefig(fig_path + subID + '_ecg_topo_eeg')

      answer = raw_input("Do you want to exclude ECG component(s) at this time? (y/n)\n")

      if answer == 'y': # yes, we want to go ahead!
        # choose which components to mark as heartbeat artefacts
        ecg_reject = raw_input('Type the indices for the heartbeat component(s) you want to reject separated by commas (e.g. 0, 3):\n')
        ecg_reject = map(int,ecg_reject.split(',')) # convert string input to a list of integers
        # add choice(s) to exclusion list using .extend
        ica.exclude.extend(ecg_reject)
      else:
        print "OK, not excluding any ECG components."
        ecg_reject = None

    else:
      print "No ECG components found! Maybe you should try tweaking the parameters."
      ecg_reject = None

  else:
    print "###\nNOT doing ECG.\n###"

  ################# Reject any additional components ###############

  # If the algorithm failed to detect either 
  if (eog_reject is None) or (ecg_reject is None) or (ica.exclude is None):

    answer = raw_input("Do you want to exclude any ECG or EOG components that weren't detected? (y/n)\n")

    if answer == 'y': # yes, we want to go ahead!
      # choose which components to mark as eyeblink artefacts
      comp_reject = raw_input('\nType the indices for the component(s) you want to reject separated by commas (e.g. 0, 3):\n')
      comp_reject = map(int,comp_reject.split(',')) # convert string input to a list of integers
      # put choice(s) in exclusion list 
      ica.exclude.extend(comp_reject)
    else:
        print "OK, not excluding any additional components."
  
  ####################### Look at the result #######################

  answer = raw_input("Have you found any components that you want to remove for this dataset? (y/n)\n")

  if answer == 'y':
    print "OK, moving on. You can still choose to not save the ICAed data at the end if you change your mind.\n"

  # estimate average artifact
  # using eog and ecg epochs created earlier and not used (could have used epochs for the artifact detection, but used raw instead)

    if do_eog =='yes' and eog_reject is not None:
      eog_evoked = eog_epochs.average()
      fig = ica.plot_sources(eog_evoked, exclude=eog_reject)  # plot EOG sources + selection
      plt.savefig(fig_path + subID + '_eog_reconstruct')


      fig = ica.plot_overlay(eog_evoked, exclude=eog_reject)  # plot EOG cleaning
      plt.savefig(fig_path + subID + '_eog_cleaning')

    if do_ecg == 'yes' and ecg_reject is not None:
      ecg_evoked = ecg_epochs.average()
      fig = ica.plot_sources(ecg_evoked, exclude=ecg_reject)  # plot ECG sources + selection
      plt.savefig(fig_path + subID + '_ecg_reconstruct')

      fig = ica.plot_overlay(ecg_evoked, exclude=ecg_reject)  # plot ECG cleaning
      plt.savefig(fig_path + subID + '_ecg_cleaning')

    # check the amplitudes do not change
    ica.plot_overlay(raw)  # EOG artifacts remain
    plt.savefig(fig_path + subID + '_ICA_clean')

    plt.show() # show plots

  ####################### Save and apply ICA solution #######################

    answer_2 = raw_input("Are you sure you want to exclude the component(s) and save? (y/n)\n")

    if answer_2 == 'y': # yes, we want to go ahead!
      raw = Raw(fname_raw_in, preload=True)

      print "###\nOK, applying ICA to raw data and rejecting component(s) {}...\n###" .format(ica.exclude)
      raw_ica = ica.apply(raw) # apply ICA to raw data

      print "###\nSaving ICAed data to " + fname_raw_out + "\n###"
      raw_ica.save(fname_raw_out, overwrite=True) # save ICAed raw data

    elif answer_2 == 'n': # no, OK - let's find out if you want to type them in again or just not save anything

      answer_3 = raw_input("OK, would you like to type in the EOG and ECG components you want to exclude again, apply the solution and save? (y/n)\n")

      if answer_3 == 'y':

        ica.exclude = [] # remove previously specified components

        raw = Raw(fname_raw_in, preload=True)

        comp_reject_2 = raw_input("Type in the component(s) you would like to reject separated by commas (e.g. 0, 2):\n")
        comp_reject_2 = map(int,comp_reject_2.split(',')) # convert string input to a list of integers

        # add choice(s) to exclusion list
        ica.exclude.extend(comp_reject_2)

        print "###\nOK, applying ICA to raw data and rejecting components {}...###\n" .format(ica.exclude)
        raw_ica = ica.apply(raw) # apply ICA to raw data

        print "###\nSaving ICAed data to " + fname_raw_out + "\n###"
        raw_ica.save(fname_raw_out, overwrite=True) # save ICAed raw data

      elif answer_3 == 'n':
        print("OK, no problem, we'll leave it for now :)")

  # if no components to reject
  elif answer == 'n':
    print "###\nNo components rejected for subject " + ss[0] + ".\n###"

  # shall we close the figures before moving on to the next subject?
  answer_4 = raw_input("Close all figures? (y/n)\n")
  if answer_4 == 'y':
    plt.close('all')
  elif answer_4 == 'n':
    print "OK, keeping figures."

  print "###\nDone with subject " + ss[0] + ".\n###\n\n"