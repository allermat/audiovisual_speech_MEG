# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:03:35 2020

@author: ma09
"""
import os
import os.path as op
import socket
import mne
import numpy as np
from mne.time_frequency import csd_morlet, csd_fourier
from mne.beamformer import make_dics, apply_dics_csd
from argparse import ArgumentParser
from scipy.io import loadmat
from megcoherence_utils import random_derangement, project_dir
# import time

# Parse input
parser = ArgumentParser()
parser.add_argument('cond', help='Conditions as hyphen separated list')
parser.add_argument('-ct', '--channel_type', help='Channel type (mag or grad)',
                    default='mag')
parser.add_argument('-ca', '--csd_alg', help='Algorithm to compute CSD (morlet or fourier)',
                    default='fourier')
parser.add_argument('-ms', '--match_sample', help='Subsample the trials to match other condition',
                    default='no')
parser.add_argument('-osd', '--out_subdir', help='Subdirectory to dave data',
                    default='')
parser.add_argument('-pc', '--part_coh',
                    help='Compute partial coherence (yes or no)', default='no')
parser.add_argument('-sc', '--square_coh', help='Compute magnitude squared coherence (yes or no)',
                    default='no')
parser.add_argument('-ss', '--source_space', help='Source space type (surf or vol)',
                    default='surf')

args = parser.parse_args()
cond = args.cond
csd_alg = args.csd_alg
channel_type = args.channel_type
part_coh = args.part_coh
out_subdir = args.out_subdir
square_coh = args.square_coh
source_space = args.source_space
match_sample = args.match_sample

# Unpack and check input
temp = cond.split('-')
if len(temp) == 2:
    mod_str, cond_str = temp
    perm_str = 'true'
elif len(temp) == 3:
    mod_str, cond_str, perm_str = temp
else:
    raise ValueError()

if cond_str == 'allAud':
    cond_sel = ['AO20', 'AO70', 'AV20', 'AV70']
elif cond_str == 'allVis':
    cond_sel = ['AV20', 'AV70', 'VO']
elif cond_str == 'AO':
    cond_sel = ['AO20', 'AO70']
elif cond_str == 'AV':
    cond_sel = ['AV20', 'AV70']
elif cond_str == 'Ahigh':
    cond_sel = ['AO70', 'AV70']
elif cond_str == 'Alow':
    cond_sel = ['AO20', 'AV20']
elif any(cond_str == x for x in ['AO20', 'AO70', 'AV20', 'AV70', 'VO']):
    cond_sel = [cond_str]
else:
    raise ValueError()

assert any(perm_str == x for x in ['true', 'perm'])

assert any(part_coh == x for x in ['yes', 'no'])

assert any(square_coh == x for x in ['yes', 'no'])

assert any(csd_alg == x for x in ['morlet', 'fourier'])

# Set which channel to partialize if required. This is always the other
# envelope channel than the one used for computing coherence.
if part_coh == 'yes':
    coh_str = 'pcoh'
    if mod_str == 'aud':
        part_mod_str = 'vis'
    elif mod_str == 'vis':
        part_mod_str = 'aud'
else:
    coh_str = 'coh'
    part_mod_str = 'none'

# Frequencies at which the coherence analysis is performed
freq = list(range(2, 7))
freq_str = '{}-{}Hz'.format(freq[0], freq[-1])

# Set how many jobs csd calculation should use
if 'node' in socket.gethostname():
    N_JOBS_CSD = 16
else:
    N_JOBS_CSD = 8

subIDlist = ['sub-' + str(n).zfill(2) for n in list(range(1, 15))]

subjects_dir = op.join(project_dir, 'data', 'derivatives', 'anat')
# # load fsaverage source space for morphing
# mne.datasets.fetch_fsaverage(subjects_dir)  # ensure fsaverage src exists
# fname_fs_src = subjects_dir + '/fsaverage/bem/fsaverage-vol-5-src.fif'
# src_fs = mne.read_source_spaces(fname_fs_src)


##############################################################################
# Define functions
def load_epochs_env(subID, dest_dir):
    ft_epochs_fname = op.join(dest_dir, 'trials_no_env_MEG.mat')
    envelopes_fname = op.join(dest_dir, 'envelopes.mat')
    raw_data_folder = op.join(project_dir, 'derivatives', 'maxfilter', subID)

    raw_file = op.join(raw_data_folder, 'concatenated_icaed_raw_trans.fif')
    raw = mne.io.read_raw_fif(raw_file, verbose=False)

    epochs = mne.read_epochs_fieldtrip(ft_epochs_fname, raw.info,
                                       data_name='trials_no_env_MEG')
    event_dict = {'AO20': 16, 'AO70': 32, 'AV20': 48, 'AV70': 64, 'VO': 80}
    epochs.event_id = event_dict

    # Load envelopes and add them to the epoched data as 'misc' channels
    temp = loadmat(envelopes_fname)
    envelopes = temp['envelopes']
    info_scratch = mne.create_info(['audEnv', 'visEnv'], 250.0,
                                   ch_types=['misc', 'misc'])
    info_env = epochs.info.copy()
    info_env['ch_names'] += info_scratch['ch_names']
    info_env['chs'] += info_scratch['chs']
    info_env['nchan'] += info_scratch['nchan']

    epochs_data = epochs.get_data()
    epochs_data = np.concatenate((epochs_data, envelopes), axis=1)
    epochs_env = mne.EpochsArray(data=epochs_data, info=info_env,
                                 events=epochs.events,
                                 event_id=epochs.event_id)
    epochs_env.apply_baseline((None, 1))

    return epochs_env


def partialize_csd(csd, part_mod_str):
    # Partialize the cross-spectral density matrix
    # Input:
    #     csd: instance of csd object
    #     part_mod_str: name of the channel to be partialized
    # Details:
    #     Based on ft_connectivity_corr.m from FieldTrip by Jan-Mathijs
    #     Schoffelen
    #     Partial spectra are computed as in Rosenberg JR et al (1998)
    #     J.Neuroscience Methods, equation 38
    csd_part = csd.copy()
    if part_mod_str == 'aud':
        penv_idx = csd.ch_names.index('audEnv')
    else:
        penv_idx = csd.ch_names.index('visEnv')
    n_chan = len(csd.ch_names)
    n_freq = len(csd.frequencies)
    # pre-allocating new data array
    data_new = np.empty((int((n_chan**2 + n_chan) / 2), n_freq))
    data_new[:] = np.nan
    for i, f in enumerate(csd.frequencies):
        # Partializing separately each frequency
        temp = csd.get_data(f)
        AA = temp
        AB = temp[:, penv_idx, None]
        BA = temp[penv_idx, :, None].transpose()
        BB = temp[penv_idx, penv_idx]
        A = AA - np.dot((AB / BB), BA)
        # Converting the partialized csd matrix to upper triangular
        # vector format
        data_new[:, i] = A[np.triu_indices_from(A)]
    csd_part._data = data_new

    return csd_part


def sensor_coherence_env(epochs, csd_meg_env, mod_str, channel_type, square_coh):
    # estimate coherence between the meg sensors
    # and the external sensor. The equation for coherence is:
    #
    #                                             |MEG-EXTERNAL CSD|
    #          MEG-EXTERNAL coherence = --------------------------------------
    #                                   sqrt(MEG POWER) * sqrt(EXTERNAL POWER)
    #
    # We can use the :func:`csd_morlet` function compute the full
    # sensor-to-sensor CSD, which contains what we need to compute
    # the numerator of the equation. The diagonal of this CSD matrix
    # contains the power for each sensor, which we can use to compute
    # the denominator of the equation.
    # To keep it consistent with FieldTrip I compute the absolute value of
    # coherency instead of the magnitude-squared coherence. Later it can be
    # squared. 

    csd_data = csd_meg_env.mean().get_data()
    # Sensor-level coherence
    psd = np.diag(csd_data).real
    coh_mat = np.abs(csd_data) / (np.sqrt(psd[np.newaxis, :]) * np.sqrt(psd[:, np.newaxis]))
    if square_coh == 'yes':
        coh_mat = np.square(coh_mat)
    freq = csd_meg_env.frequencies
    # save topomap of coherence
    if mod_str == 'aud':
        env_idx = csd_meg_env.ch_names.index('audEnv')
    else:
        env_idx = csd_meg_env.ch_names.index('visEnv')
    coh_env = coh_mat[:-2, env_idx]
    info_coh = mne.pick_info(epochs_env.info,
                             mne.pick_types(epochs.info, meg=channel_type))
    coh = mne.time_frequency.AverageTFR(info_coh,
                                        coh_env[:, np.newaxis, np.newaxis],
                                        np.array([0]),
                                        np.array([np.mean(freq)]),
                                        len(freq))
    
    return coh


def dics_coherence_env(csd, dics, mod_str, square_coh):
    # To estimate source-level coherence with the external sensor, we can use the
    # DICS beamformer to::
    #  1. Project the CSD between each gradiometer and the external source to the
    #     cortical surface to compute the numerator of the equation.
    #  2. Project the gradiometer CSD to the cortical surface to compute the
    #     denominator of the equation.
    #
    if mod_str == 'aud':
        env_idx = csd.ch_names.index('audEnv')
    else:
        env_idx = csd.ch_names.index('visEnv')

    source_power, freq = apply_dics_csd(csd, dics, verbose=False)
    coherence = []
    for freq_idx, act_freq in enumerate(freq):
        # Compute source coherence separately for each frequency
        act_source_power = source_power.data[:, freq_idx]
        csd_data = csd.get_data(act_freq)
        psd = np.diag(csd_data).real
        external_power = psd[env_idx]

        source_csd = dics['weights'][freq_idx].dot(csd_data[:-2, env_idx])
        act_coh = np.absolute(source_csd) / (np.sqrt(act_source_power) * np.sqrt(external_power))
        if square_coh == 'yes':
            act_coh = np.square(act_coh)
        coherence.append(act_coh[np.newaxis, :])
    coherence = np.concatenate(coherence, axis=0).mean(axis=0)
    # Average source coherence across frequencies
    stc_coh = mne.SourceEstimate(coherence[:, np.newaxis],
                                 vertices=dics['vertices'], tmin=0, tstep=1)
    return stc_coh


def permute_envelopes(epochs):
    epochs_data = epochs.get_data()
    envelopes = epochs_data[:, [-2, -1], :]
    # Permute the envelopes across trials using random derangements
    idx = random_derangement(envelopes.shape[0])
    envelopes = envelopes[idx, :, :]
    epochs_data[:, [-2, -1], :] = envelopes
    epochs._data = epochs_data

    return epochs


def subsample_epochs(epochs, n_resample):
    # Making sure that all conditions are roughly equally represented
    # in the subsampled epochs array
    event_ids = epochs.events[:, 2]
    unique_events, event_counts = np.unique(event_ids, return_counts=True)
    n_cond = len(unique_events)
    event_counts_new = np.floor(event_counts/n_cond).astype(int)
    if event_counts_new.sum() != n_resample:
        # In case event counts were not divisible with n_cond,
        # the number of selected trials will not match n_resample, but
        # will be off by 1. I add 1 trial to a randomly chosen condition,
        # this way over many repetitions the average number of trials after
        # subsampling will be roughly uniform across conditions.
        event_counts_new[np.random.choice(n_cond)] += 1
    # Subsample trials separately in each condition
    idx_list = []
    for i, cond_id in enumerate(unique_events):
        act_idx = np.nonzero(event_ids == cond_id)[0]
        idx_list.append(np.random.choice(act_idx, size=event_counts_new[i],
                                         replace=False))
    # Drop the selected epochs
    idx = np.concatenate(idx_list)
    epochs = epochs.drop(idx)

    return epochs


def process_subject(subID, epochs, forward, mod_str, freq, channel_type,
                    perm_str, part_mod_str, csd_alg, square_coh, match_sample):
    forward = mne.pick_types_forward(forward, meg=channel_type, eeg=False)

    if match_sample == 'no':  # Don't subsample epochs
        # Estimate the cross-spectral density (CSD) matrix
        # t = time.time()
        if csd_alg == 'morlet':
            csd_meg = csd_morlet(epochs[cond_sel], frequencies=freq,
                                 picks=[channel_type], verbose=False,
                                 n_jobs=N_JOBS_CSD)
        elif csd_alg == 'fourier':
            # To make sure min and max frequencies are included I use a
            # little offset
            offset = 0.2
            csd_meg = csd_fourier(epochs[cond_sel], fmin=freq[0]-offset,
                                  fmax=freq[-1]+offset,
                                  picks=[channel_type], verbose=False,
                                  n_jobs=N_JOBS_CSD)
        # print('Elapsed time: {}'.format(time.time() - t))
        # Compute the DICS power map.
        dics = make_dics(epochs.info, forward, csd_meg, reg=0.05,
                         pick_ori='max-power')
        if perm_str == 'true':  # compute true coherence
            # Estimate the CSD matrix with envelopes included
            if csd_alg == 'morlet':
                csd_meg_env = csd_morlet(epochs[cond_sel], frequencies=freq,
                                         picks=[channel_type, 'misc'], verbose=False,
                                         n_jobs=N_JOBS_CSD)
            elif csd_alg == 'fourier':
                csd_meg_env = csd_fourier(epochs[cond_sel], fmin=freq[0]-offset,
                                          fmax=freq[-1]+offset,
                                          picks=[channel_type, 'misc'],
                                          verbose=False, n_jobs=N_JOBS_CSD)
            # Partialize the csd matrix if partial coherence is required
            if part_mod_str != 'none':
                csd_meg_env = partialize_csd(csd_meg_env, part_mod_str)

            # Sensor-level coherence with external source
            coh = sensor_coherence_env(epochs, csd_meg_env, mod_str,
                                       channel_type, square_coh)
            # Source-level coherence with external source
            stc_coh = dics_coherence_env(csd_meg_env, dics, mod_str, square_coh)
        else:  # compute permuted coherence
            # Setup the random number generator
            rng = np.random.RandomState()
            n_perm = 100
            stc_data_list, coh_data_list = [], []
            for i in range(0, n_perm):
                epochs_perm = permute_envelopes(epochs[cond_sel])
                # Estimate the CSD matrix with envelopes included
                if csd_alg == 'morlet':
                    csd_meg_env = csd_morlet(epochs_perm, frequencies=freq,
                                             picks=[channel_type, 'misc'],
                                             verbose=False, n_jobs=N_JOBS_CSD)
                elif csd_alg == 'fourier':
                    csd_meg_env = csd_fourier(epochs_perm, fmin=freq[0]-offset,
                                              fmax=freq[-1]+offset,
                                              picks=[channel_type, 'misc'],
                                              verbose=False, n_jobs=N_JOBS_CSD)
                # Partialize the csd matrix if partial coherence is required
                if part_mod_str != 'none':
                    csd_meg_env = partialize_csd(csd_meg_env, part_mod_str)
                # Sensor-level coherence with external source
                coh_temp = sensor_coherence_env(epochs_perm, csd_meg_env, mod_str,
                                                channel_type, square_coh)
                coh_data_list.append(coh_temp.data)
                # Source-level coherence with external source
                stc_coh_temp = dics_coherence_env(csd_meg_env, dics, mod_str,
                                                  square_coh)
                stc_data_list.append(stc_coh_temp.data)

            coh = coh_temp.copy()
            coh._data = np.concatenate(coh_data_list, axis=2).mean(
                axis=2, keepdims=True)
            stc_coh = stc_coh_temp.copy()
            stc_coh._data = np.concatenate(stc_data_list, axis=1).mean(
                axis=1, keepdims=True)

    else:  # subsample epochs
        # Check to how many trials must the data be subsampled
        n_subsample = len(epochs[match_sample])
        # Setup the random number generator
        rng = np.random.RandomState()
        n_rep = 100
        stc_data_list, coh_data_list = [], []
        for i in range(0, n_rep):
            # If resampling is done the DICS object must also be computed
            # on the resampled data
            # Estimate the cross-spectral density (CSD) matrix
            act_epochs = subsample_epochs(epochs[cond_sel], n_subsample)
            if csd_alg == 'morlet':
                csd_meg = csd_morlet(act_epochs, frequencies=freq,
                                     picks=[channel_type], verbose=False,
                                     n_jobs=N_JOBS_CSD)
            elif csd_alg == 'fourier':
                # To make sure min and max frequencies are included I use a
                # little offset
                offset = 0.2
                csd_meg = csd_fourier(act_epochs, fmin=freq[0]-offset,
                                      fmax=freq[-1]+offset,
                                      picks=[channel_type], verbose=False,
                                      n_jobs=N_JOBS_CSD)
                # Compute the DICS power map.
            dics = make_dics(act_epochs.info, forward, csd_meg, reg=0.05,
                             pick_ori='max-power')

            if perm_str == 'perm':  # shuffle envelopes if required
                act_epochs = permute_envelopes(act_epochs)
            # Estimate the CSD matrix with envelopes included
            if csd_alg == 'morlet':
                csd_meg_env = csd_morlet(act_epochs, frequencies=freq,
                                         picks=[channel_type, 'misc'],
                                         verbose=False, n_jobs=N_JOBS_CSD)
            elif csd_alg == 'fourier':
                csd_meg_env = csd_fourier(act_epochs, fmin=freq[0]-offset,
                                          fmax=freq[-1]+offset,
                                          picks=[channel_type, 'misc'],
                                          verbose=False, n_jobs=N_JOBS_CSD)
            # Partialize the csd matrix if partial coherence is required
            if part_mod_str != 'none':
                csd_meg_env = partialize_csd(csd_meg_env, part_mod_str)
            # Sensor-level coherence with external source
            coh_temp = sensor_coherence_env(act_epochs, csd_meg_env, mod_str,
                                            channel_type, square_coh)
            coh_data_list.append(coh_temp.data)
            # Source-level coherence with external source
            stc_coh_temp = dics_coherence_env(csd_meg_env, dics, mod_str,
                                              square_coh)
            stc_data_list.append(stc_coh_temp.data)

        coh = coh_temp.copy()
        coh._data = np.concatenate(coh_data_list, axis=2).mean(
            axis=2, keepdims=True)
        stc_coh = stc_coh_temp.copy()
        stc_coh._data = np.concatenate(stc_data_list, axis=1).mean(
            axis=1, keepdims=True)

    return coh, stc_coh


##############################################################################
# Execute across subjects
for subID in subIDlist:
    source_dir = op.join('project_dir, 'data', 'derivatives', 'megcoherence', subID)
    if out_subdir:
        dest_dir = op.join(source_dir, out_subdir)
        if not op.exists(dest_dir):
            os.mkdir(dest_dir)
    else:
        dest_dir = source_dir
    
    # Load epochs with envelopes
    epochs_env = load_epochs_env(subID, source_dir)

    # Load precomputed forward solution
    forward = mne.read_forward_solution(op.join(source_dir, subID+'_surf-fwd.fif'))
    forward = mne.convert_forward_solution(forward, surf_ori=True)

    # Compute sensor and source level coherence with the envelope
    coh, stc_coh = process_subject(subID, epochs_env, forward, mod_str, freq,
                                   channel_type, perm_str, part_mod_str,
                                   csd_alg, square_coh, match_sample)

    # Save coherence estimates
    if match_sample != 'no':
        match_str = '_match{}'.format(match_sample)
    else:
        match_str = ''
    coh_fname = '{}_mne_{}_{}-{}_{}_{}_{}{}-tfr.h5'.format(
        subID, coh_str, mod_str, cond_str, freq_str, channel_type, perm_str,
        match_str)
    coh.save(op.join(dest_dir, coh_fname), overwrite=True)
    stc_fname = '{}_mne_src-{}_{}-{}_{}_{}_{}{}'.format(
        subID, coh_str, mod_str, cond_str, freq_str, channel_type, perm_str,
        match_str)
    stc_coh.save(op.join(dest_dir, stc_fname))
