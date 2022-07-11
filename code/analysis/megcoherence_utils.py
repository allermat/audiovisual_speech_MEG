# -*- coding: utf-8 -*-
"""
Utility functions for megcoherence analysis

Created on Wed Dec 02 16:19:35 2020

@author: Máté Aller, allermat@gmail.com
"""

import os.path as op
from os import listdir
import mne
import numpy as np
import pandas as pd
import re


def load_data_roi_analysis(conds, rois, channel_type='mag', freqs=[2, 6], save=False, sub_dir=None):
    """ Load data for anatomical ROI analysis to a pandas dataframe"""

    subjects_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                           'data', 'derivatives', 'anat')
    if sub_dir:
        dest_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                           'data', 'derivatives', 'megcoherence', 'group',
                           sub_dir)
    else:
        dest_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                           'data', 'derivatives', 'megcoherence', 'group')
    subj_dir_list = [f for f in listdir(subjects_dir)
                     if op.isdir(op.join(subjects_dir, f))]
    subj_id_list = sorted([s for s in subj_dir_list if s.startswith('sub')])
    n_subjects = len(subj_id_list)
    freq_str = '{}-{}Hz'.format(freqs[0], freqs[1])

    # Load anatomical labels for rois
    anat_labels = load_label_roi(rois, subjects_dir)
    # Initialize lists for collecting data
    sub, acc_clar, stim_mod, coh_mod, part, perm, match, roi, hemi, coh_mean, \
        coh_max = ([] for i in range(11))
    for act_cond in conds:

        cond_str, cond_mod, cond_perm, cond_coh_str, cond_stim_mod, \
            cond_match, cond_acc_clar = process_cond(act_cond)
        temp_mean, temp_max = [], []
        for subID in subj_id_list:
            if sub_dir:
                source_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                                     'data', 'derivatives', 'megcoherence', subID, sub_dir)
            else:
                source_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                                     'data', 'derivatives', 'megcoherence', subID)
            evoked_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                                 'data', 'derivatives', 'megevoked', subID)
            # Load precomputed forward solution and source estimates
            forward_fname = op.join(evoked_dir, ''.join([subID, '_surf-fwd.fif']))
            forward = mne.read_forward_solution(forward_fname)

            if cond_match:
                temp_cond_match = '_'+cond_match
            else:
                temp_cond_match = cond_match
            stc_fname = '{}_mne_src-{}_{}-{}_{}_{}_{}{}-lh.stc'.format(
                subID, cond_coh_str, cond_mod, cond_str, freq_str,
                channel_type, cond_perm, temp_cond_match)
            stc = mne.read_source_estimate(op.join(source_dir, stc_fname))

            # Morph source estimates to fsaverage
            morph = mne.compute_source_morph(
                forward['src'], subject_from=subID, subject_to='fsaverage',
                subjects_dir=subjects_dir)
            stc_fs = morph.apply(stc)

            stc_data = []
            for lab in anat_labels:
                # Apply label
                stc_fs_roi = stc_fs.in_label(lab)
                # Collect source estimates
                stc_data.append(stc_fs_roi.data)
            # Compute mean across vertices in each ROI
            temp = [np.mean(x, keepdims=True) for x in stc_data]
            # Concatenate the means across rois, each element of the list
            # temp_mean is 1 x nROIs
            temp_mean.append(np.concatenate(temp, axis=1))
            temp = [np.amax(x, keepdims=True) for x in stc_data]
            temp_max.append(np.concatenate(temp, axis=1))

        n_rep = n_subjects*len(anat_labels)
        sub.append(np.tile(subj_id_list, len(anat_labels)))
        coh_mean.append(np.concatenate(temp_mean, axis=0).flatten('F'))
        coh_max.append(np.concatenate(temp_max, axis=0).flatten('F'))
        acc_clar.append(np.repeat(cond_acc_clar, n_rep))
        stim_mod.append(np.repeat(cond_stim_mod, n_rep))
        coh_mod.append(np.repeat(cond_mod, n_rep))
        part.append(np.repeat(cond_coh_str == 'pcoh', n_rep))
        perm.append(np.repeat(cond_perm, n_rep))
        if not cond_match:
            match.append(np.repeat('none', n_rep))
        else:
            match.append(np.repeat(cond_match, n_rep))
        roi.append(np.repeat(rois, n_subjects))
        hemi.append(np.repeat('both', n_rep))

    # Arranging data in a dataframe
    df = pd.DataFrame(
        {'subID': pd.Categorical(np.concatenate(sub)),
         'partialized': pd.Categorical(np.concatenate(part)),
         'coh_modality': pd.Categorical(np.concatenate(coh_mod)),
         'stim_modality': pd.Categorical(np.concatenate(stim_mod)),
         'acc_clarity': pd.Categorical(np.concatenate(acc_clar)),
         'permutation': pd.Categorical(np.concatenate(perm)),
         'match_sample': pd.Categorical(np.concatenate(match)),
         'roi': pd.Categorical(np.concatenate(roi)),
         'hemisphere': pd.Categorical(np.concatenate(hemi)),
         'coh_mean': pd.Categorical(np.concatenate(coh_mean)),
         'coh_max': pd.Categorical(np.concatenate(coh_max))})

    if save:
        df_name = 'megcoherence_anatomical_roi_data.csv'
        df.to_csv(op.join(dest_dir, df_name), index=False)

    return df


def load_data_cluster(conds, rois, channel_type='mag', freqs=[2, 6], save=False, sub_dir=None):
    """ Load data from significant clusters to a pandas dataframe"""

    assert len(conds) == len(rois), 'Labels must be specified separately for each condition'

    subjects_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                           'data', 'derivatives', 'anat')
    if sub_dir:
        dest_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                           'data', 'derivatives', 'megcoherence', 'group',
                           sub_dir)
    else:
        dest_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                           'data', 'derivatives', 'megcoherence', 'group')
    subj_dir_list = [f for f in listdir(subjects_dir)
                     if op.isdir(op.join(subjects_dir, f))]
    subj_id_list = sorted([s for s in subj_dir_list if s.startswith('sub')])
    n_subjects = len(subj_id_list)
    freq_str = '{}-{}Hz'.format(freqs[0], freqs[1])

    # Load anatomical labels for rois
    func_labels = load_label_func(rois, dest_dir, subjects_dir)

    # Initialize lists for collecting data
    sub, acc_clar, stim_mod, coh_mod, part, perm, match, roi, hemi, coh_mean, \
        coh_max = ([] for i in range(11))
    for i, act_cond in enumerate(conds):

        cond_str, cond_mod, cond_perm, cond_coh_str, cond_stim_mod, \
            cond_match, cond_acc_clar = process_cond(act_cond)
        temp_mean, temp_max, temp_hemi = [], [], []
        for subID in subj_id_list:
            if sub_dir:
                source_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                                     'data', 'derivatives', 'megcoherence', subID, sub_dir)
            else:
                source_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                                     'data', 'derivatives', 'megcoherence', subID)
            evoked_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                                 'data', 'derivatives', 'megevoked', subID)
            # Load precomputed forward solution and source estimates
            forward_fname = op.join(evoked_dir, ''.join([subID, '_surf-fwd.fif']))
            forward = mne.read_forward_solution(forward_fname)

            if cond_match:
                temp_cond_match = '_'+cond_match
            else:
                temp_cond_match = cond_match
            stc_fname = '{}_mne_src-{}_{}-{}_{}_{}_{}{}-lh.stc'.format(
                subID, cond_coh_str, cond_mod, cond_str, freq_str,
                channel_type, cond_perm, temp_cond_match)
            stc = mne.read_source_estimate(op.join(source_dir, stc_fname))

            # Morph source estimates to fsaverage
            morph = mne.compute_source_morph(
                forward['src'], subject_from=subID, subject_to='fsaverage',
                subjects_dir=subjects_dir)
            stc_fs = morph.apply(stc)

            stc_data = []
            for lab in func_labels[i]:
                # Apply label
                stc_fs_roi = stc_fs.in_label(lab)
                # Collect source estimates
                stc_data.append(stc_fs_roi.data)
            # Compute mean across vertices in each ROI
            temp = [np.mean(x, keepdims=True) for x in stc_data]
            # Concatenate the means across rois, each element of the list
            # temp_mean is 1 x nROIs
            temp_mean.append(np.concatenate(temp, axis=1))
            temp = [np.amax(x, keepdims=True) for x in stc_data]
            temp_max.append(np.concatenate(temp, axis=1))
            temp = [x.hemi for x in func_labels[i]]
            temp_hemi.append(np.expand_dims(temp, axis=0))
            # temp_hemi.append([temp])

        n_rep = n_subjects*len(func_labels[i])
        sub.append(np.tile(subj_id_list, len(func_labels[i])))
        coh_mean.append(np.concatenate(temp_mean, axis=0).flatten('F'))
        coh_max.append(np.concatenate(temp_max, axis=0).flatten('F'))
        acc_clar.append(np.repeat(cond_acc_clar, n_rep))
        stim_mod.append(np.repeat(cond_stim_mod, n_rep))
        coh_mod.append(np.repeat(cond_mod, n_rep))
        part.append(np.repeat(cond_coh_str == 'pcoh', n_rep))
        perm.append(np.repeat(cond_perm, n_rep))
        if not cond_match:
            match.append(np.repeat('none', n_rep))
        else:
            match.append(np.repeat(cond_match, n_rep))
        roi.append(np.repeat(rois[i], n_rep))
        hemi.append(np.concatenate(temp_hemi, axis=0).flatten('F'))

    # Arranging data in a dataframe
    df = pd.DataFrame(
        {'subID': pd.Categorical(np.concatenate(sub)),
         'partialized': pd.Categorical(np.concatenate(part)),
         'coh_modality': pd.Categorical(np.concatenate(coh_mod)),
         'stim_modality': pd.Categorical(np.concatenate(stim_mod)),
         'acc_clarity': pd.Categorical(np.concatenate(acc_clar)),
         'permutation': pd.Categorical(np.concatenate(perm)),
         'match_sample': pd.Categorical(np.concatenate(match)),
         'roi': pd.Categorical(np.concatenate(roi)),
         'hemisphere': pd.Categorical(np.concatenate(hemi)),
         'coh_mean': pd.Categorical(np.concatenate(coh_mean)),
         'coh_max': pd.Categorical(np.concatenate(coh_max))})

    if save:
        df_name = 'megcoherence_functional_roi_data.csv'
        df.to_csv(op.join(dest_dir, df_name), index=False)

    return df


def process_cond(cond):
    # Check input and set variable accordingly
    temp = cond.split('-')
    if len(temp) == 2:
        coh_mod, cond_str = temp
        perm_str = 'true'
        coh_str = 'coh'
        match_str = ''
    elif len(temp) == 3:
        coh_mod, cond_str, perm_str = temp
        coh_str = 'coh'
        match_str = ''
    elif len(temp) == 4:
        coh_str, coh_mod, cond_str, perm_str = temp
        match_str = ''
    elif len(temp) == 5:
        coh_str, coh_mod, cond_str, perm_str, match_str = temp
    else:
        raise ValueError()
    assert any(coh_str == x for x in ['coh', 'pcoh']), 'Must be coh or pcoh!'

    # Further process cond_str
    p = re.compile('([A-Z]{2})([0-9]{2})?')
    m = p.match(cond_str)
    if m:
        # This will match all but Alow, Ahigh
        stim_mod = m.group(1)
        if m.group(2) == '20':
            acc_clar = 'low'
        elif m.group(2) == '70':
            acc_clar = 'high'
        elif m.group(2) is None:
            acc_clar = 'none'
        else:
            raise ValueError()
    else:
        if any(cond_str == s for s in ['Alow', 'Ahigh']):
            stim_mod = 'A'
            acc_clar = cond_str.split('A')[-1]
        elif any(cond_str == s for s in ['allAud', 'allVis']):
            stim_mod = cond_str[0:4]
            acc_clar = 'none'
        else:
            raise ValueError()

    return cond_str, coh_mod, perm_str, coh_str, stim_mod, match_str, acc_clar


def load_label_roi(rois, subjects_dir, hemi='both'):
    # Load anatomical labels for ROIs
    label = []
    for roi in rois:
        if roi == 'Occ_cortex':
            anat_label = [
                'G_and_S_occipital_inf', 'G_cuneus', 'G_occipital_middle',
                'G_occipital_sup', 'Pole_occipital',
                'S_oc_middle_and_Lunatus', 'S_oc_sup_and_transversal']
        elif roi == 'STG':
            anat_label = ['G_temp_sup-Lateral']
        else:
            anat_label = [roi]

        # Load label
        def fun_read_label(lab, h):
            return mne.read_labels_from_annot('fsaverage', parc='aparc.a2009s',
                                              hemi=h, subjects_dir=subjects_dir,
                                              regexp=lab)[0]

        if hemi == 'both':
            temp = [fun_read_label(s, 'lh') for s in anat_label]
            temp_label_lh = temp[0]
            for x in temp[1:]:
                temp_label_lh += x
            temp = [fun_read_label(s, 'rh') for s in anat_label]
            temp_label_rh = temp[0]
            for x in temp[1:]:
                temp_label_rh += x
            label.append(temp_label_lh + temp_label_rh)
        else:
            temp = [fun_read_label(s, hemi) for s in anat_label]
            label.append(temp[0])
            for x in temp[1:]:
                label[-1] += x

    return label


def load_label_func(func_rois, source_dir, subjects_dir):
    # Load functinally defined ROIs as labels
    label = []
    for act_roi in func_rois:

        if act_roi == 'vis_VO>AV':
            func_label_fname = 'group_src-pcoh-clust_vis-VO-true_-_vis-VO-perm_vs_vis-AV-true-matchVO_-_vis-AV-perm-matchVO_2-6Hz_mag-rh.stc'

        elif act_roi == 'aud_AV>AO':
            func_label_fname = 'group_src-pcoh-clust_aud-AV-true_-_aud-AV-perm_vs_aud-AO-true_-_aud-AO-perm_2-6Hz_mag-rh.stc'
        else:
            raise ValueError()

        func_label_stc = mne.read_source_estimate(op.join(source_dir, func_label_fname))
        fs_src_fname = subjects_dir + '/fsaverage/bem/fsaverage-ico-5-src.fif'
        src_fs = mne.read_source_spaces(fs_src_fname)
        func_label = mne.stc_to_label(func_label_stc, src=src_fs, smooth=True,
                                      connected=False, subjects_dir=subjects_dir)
        # Get rid of None label in case there was no cluster in a hemisphere
        func_label = [label for label in func_label if label is not None]
        label.append(func_label)
    return label


def random_derangement(n):
    """Generate a random derangement of array range
    
    This algorithm is based on https://uk.mathworks.com/matlabcentral/fileexchange/30189-randpermfull
    See also: https://stackoverflow.com/questions/25200220/generate-a-random-derangement-of-a-list
    
    Parameters
    ----------
    n : int
        Randomly permute np.arange(n) such that no elements remain at their
            original index
    
    Returns
    -------
    v : numpy ndarray
        Permutaed array range
    """
    has_violation = True
    while has_violation:
        v = np.random.permutation(n)
        has_violation = False
        for i in range(n):
            if v[i] == i:
                has_violation = True
                break
    return v


def test_random_derangement(n):
    # enumerate all derangements for testing
    import itertools
    counter = {}
    for p in itertools.permutations(range(n)):
        if all(p[i] != i for i in p):
            counter[p] = 0

    # make M probes for each derangement
    M = 5000
    for _ in range(M*len(counter)):
        # generate a random derangement
        p = tuple(random_derangement(n))
        # is it really?
        assert p in counter
        # ok, record it
        counter[p] += 1

    # the distribution looks uniform
    for p, c in sorted(counter.items()):
        print((p, c))


if __name__ == "__main__":
    df = load_data_roi_analysis(['pcoh-aud-AV-true'],
                                ['G_occipital_sup', 'Pole_occipital'])
