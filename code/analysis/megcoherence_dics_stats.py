# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:03:35 2020

@author: ma09
"""
import os.path as op
import socket
import mne
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc, permutation_t_test)
from argparse import ArgumentParser
import re
import numpy as np
from scipy import stats as stats
from scipy.io import savemat
import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.stats.anova import AnovaRM
import seaborn as sns

parser = ArgumentParser()
parser.add_argument('cond1', help='First condition')
parser.add_argument('cond2', help='Second condition')
parser.add_argument('-c3', '--cond3', help='Optinal third condition',
                    default='')
parser.add_argument('-c4', '--cond4', help='Optinal fourth condition',
                    default='')
parser.add_argument('-ct', '--channel_type', help='Channel type (mag or grad)',
                    default='mag')
parser.add_argument('-f', '--fact', help='Order of factors for rmANOVA',
                    default='')
parser.add_argument('-m', '--mask', help='Brain mask within which to do stats',
                    default='')
parser.add_argument('-he', '--hemi', help='Hemisphere to analyse (both, lh or rh)',
                    default='both')
parser.add_argument('-s', '--stat',
                    help='Statistical test type (clust, avg, max_vert)',
                    default='clust')
parser.add_argument('-sd', '--subdir', help='Subdirectory to dave data',
                    default='')
parser.add_argument('-pc', '--p_clust', help='Cluster defining p threshold',
                    default='0.05')
parser.add_argument('-ps', '--p_stat',
                    help='Statistical comparison p threshold',
                    default='0.05')
parser.add_argument('-t', '--tail',
                    help='Statistical comparison tail (-1, 0, 1)',
                    default='0')
args = parser.parse_args()
cond1 = args.cond1
cond2 = args.cond2
cond3 = args.cond3
cond4 = args.cond4
channel_type = args.channel_type
fact = args.fact.split('-')
mask = args.mask
hemi = args.hemi
stat = args.stat
subdirectory = args.subdir
p_clust = args.p_clust
p_stat = float(args.p_stat)
tail = int(args.tail)


# Check input and set variable accordingly
def process_cond(cond):
    temp = cond.split('-')
    if len(temp) == 2:
        cond_mod, cond_str = temp
        cond_perm = 'true'
        cond_coh_str = 'coh'
        cond_match_str = ''
    elif len(temp) == 3:
        cond_mod, cond_str, cond_perm = temp
        cond_coh_str = 'coh'
        cond_match_str = ''
    elif len(temp) == 4:
        cond_coh_str, cond_mod, cond_str, cond_perm = temp
        cond_match_str = ''
    elif len(temp) == 5:
        cond_coh_str, cond_mod, cond_str, cond_perm, cond_match_str = temp
    else:
        raise ValueError()
    assert any(cond_coh_str == x for x in ['coh', 'pcoh']), 'Must be coh or pcoh!'

    return cond_str, cond_mod, cond_perm, cond_coh_str, cond_match_str


def prepend_if_not_empty(str, char):
    return char+str if str else str


cond1_str, cond1_mod, cond1_perm, cond1_coh_str, cond1_match_str = process_cond(cond1)
cond2_str, cond2_mod, cond2_perm, cond2_coh_str, cond2_match_str = process_cond(cond2)
n_conds = 2
if cond3:
    assert cond4, 'cond3 and cond4 must be specified together!'
    cond3_str, cond3_mod, cond3_perm, cond3_coh_str, cond3_match_str = process_cond(cond3)
    cond4_str, cond4_mod, cond4_perm, cond4_coh_str, cond4_match_str = process_cond(cond4)
    n_conds += 2

assert any(hemi == x for x in ['both', 'lh', 'rh'])

assert any(stat == x for x in ['clust', 'avg', 'max_vert'])

freq_str = '{}-{}Hz'.format(2, 6)
views = ['lat', 'med']

# Set how many jobs csd calculation should use
if 'node' in socket.gethostname():
    N_JOBS = 16
else:
    N_JOBS = 8

n_subjects = 14

subIDlist = ['sub-' + str(n).zfill(2) for n in list(range(1, n_subjects+1))]

subjects_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects',
                       'AVSpeechMEG', 'data', 'derivatives', 'anat')
if subdirectory:
    dest_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects',
                       'AVSpeechMEG', 'data', 'derivatives', 'megcoherence',
                       'group', subdirectory)
else:
    dest_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects',
                       'AVSpeechMEG', 'data', 'derivatives', 'megcoherence',
                       'group')
anat_label = ''
# Loading mask
if not mask:
    exclude = None
else:
    p = re.compile('group_src-mask_.*')
    match = p.match(mask)
    if match:
        # Load source estimate for mask (This needs to come from a statistic
        # test created by summarize_clusters_stc() or anatomically defined ROI
        stc_mask = mne.read_source_estimate(op.join(dest_dir, mask))
        temp = stc_mask.data[:, 0]   # First cluster is the union of all
        # Select vertices which are not part of the cluster, these are going to
        # be excluded from the statistical analysis
        exclude = np.where(temp == 0)[0].tolist()
        include = np.where(temp == 1)[0].tolist()
    else:
        # The mask is anatomically defined label
        anat_label = mask
        # This only works with avg and max_vert stats
        assert any(stat == x for x in ['avg', 'max_vert']), 'Does not work with cluster statistic!'

        # Load and apply label
        if anat_label == 'Occ_cortex':
            anat_label = [
                'G_and_S_occipital_inf', 'G_cuneus', 'G_occipital_middle',
                'G_occipital_sup', 'Pole_occipital',
                'S_oc_middle_and_Lunatus', 'S_oc_sup_and_transversal']
        else:
            anat_label = [anat_label]

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
            label = temp_label_lh + temp_label_rh
        else:
            temp = [fun_read_label(s, hemi) for s in anat_label]
            label = temp[0]
            for x in temp[1:]:
                label += x

# load fsaverage source space for morphing
mne.datasets.fetch_fsaverage(subjects_dir)  # ensure fsaverage src exists
fs_src_fname = subjects_dir + '/fsaverage/bem/fsaverage-ico-5-src.fif'
src_fs = mne.read_source_spaces(fs_src_fname)
fsave_vertices = [s['vertno'] for s in src_fs]
stc_morph_list_cond1, stc_morph_list_cond2 = [], []
stc_morph_list_cond3, stc_morph_list_cond4 = [], []
for subID in subIDlist:

    if subdirectory:
        source_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                             'data', 'derivatives', 'megcoherence', subID,
                             subdirectory)
    else:
        source_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG',
                             'data', 'derivatives', 'megcoherence', subID)
    evoked_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects', 'AVSpeechMEG', 'data',
                         'derivatives', 'megevoked', subID)
    # Load precomputed forward solution and source estimates
    forward_fname = op.join(evoked_dir, ''.join([subID, '_surf-fwd.fif']))
    forward = mne.read_forward_solution(forward_fname)

    # appending an underscore to the matchstring if it is not empty
    stc_fname = '{}_mne_src-{}_{}-{}_{}_{}_{}{}-lh.stc'.format(
        subID, cond1_coh_str, cond1_mod, cond1_str, freq_str, channel_type,
        cond1_perm, prepend_if_not_empty(cond1_match_str, '_'))
    stc_cond1 = mne.read_source_estimate(op.join(source_dir, stc_fname))
    stc_fname = '{}_mne_src-{}_{}-{}_{}_{}_{}{}-lh.stc'.format(
        subID, cond2_coh_str, cond2_mod, cond2_str, freq_str, channel_type,
        cond2_perm, prepend_if_not_empty(cond2_match_str, '_'))
    stc_cond2 = mne.read_source_estimate(op.join(source_dir, stc_fname))
    if cond3:
        stc_fname = '{}_mne_src-{}_{}-{}_{}_{}_{}{}-lh.stc'.format(
            subID, cond3_coh_str, cond3_mod, cond3_str, freq_str, channel_type,
            cond3_perm, prepend_if_not_empty(cond3_match_str, '_'))
        stc_cond3 = mne.read_source_estimate(op.join(source_dir, stc_fname))
        stc_fname = '{}_mne_src-{}_{}-{}_{}_{}_{}{}-lh.stc'.format(
            subID, cond4_coh_str, cond4_mod, cond4_str, freq_str, channel_type,
            cond4_perm, prepend_if_not_empty(cond4_match_str, '_'))
        stc_cond4 = mne.read_source_estimate(op.join(source_dir, stc_fname))

    # Morph source estimates to fsaverage for group analysis
    morph = mne.compute_source_morph(
        forward['src'], subject_from=subID, subject_to='fsaverage',
        subjects_dir=subjects_dir)
    stc_cond1_fs = morph.apply(stc_cond1)
    stc_cond2_fs = morph.apply(stc_cond2)
    if cond3:
        stc_cond3_fs = morph.apply(stc_cond3)
        stc_cond4_fs = morph.apply(stc_cond4)
    if anat_label:
        # Apply label
        stc_cond1_fs = stc_cond1_fs.in_label(label)
        stc_cond2_fs = stc_cond2_fs.in_label(label)
        if cond3:
            stc_cond3_fs = stc_cond3_fs.in_label(label)
            stc_cond4_fs = stc_cond4_fs.in_label(label)
    stc_morph_list_cond1.append(stc_cond1_fs)
    stc_morph_list_cond2.append(stc_cond2_fs)
    if cond3:
        stc_morph_list_cond3.append(stc_cond3_fs)
        stc_morph_list_cond4.append(stc_cond4_fs)

# Collect source estimates across subjects into a single array
stc_morph_data_cond1 = [m.data[:, :, np.newaxis] for m in stc_morph_list_cond1]
stc_morph_data_cond2 = [m.data[:, :, np.newaxis] for m in stc_morph_list_cond2]
X_cond1 = np.concatenate(stc_morph_data_cond1, axis=2)
X_cond2 = np.concatenate(stc_morph_data_cond2, axis=2)
# Compute cond1-cond2 contrast, only for cluster statistics
X = X_cond1 - X_cond2
# If 4 conditions are specified compute interaction (cond1-cond2)-(cond3-cond4)
if cond3:
    stc_morph_data_cond3 = [m.data[:, :, np.newaxis] for m in stc_morph_list_cond3]
    stc_morph_data_cond4 = [m.data[:, :, np.newaxis] for m in stc_morph_list_cond4]
    X_cond3 = np.concatenate(stc_morph_data_cond3, axis=2)
    X_cond4 = np.concatenate(stc_morph_data_cond4, axis=2)
    X = (X_cond1 - X_cond2)-(X_cond3 - X_cond4)

# Compute statistic
# -----------------
if stat == 'clust':
    temp = p_clust.split('-')
    if len(temp) == 1:
        if tail == 0:
            thresh_clust = -stats.distributions.t.ppf(float(temp[0])/2.,
                                                      n_subjects - 1)
        else:
            thresh_clust = -stats.distributions.t.ppf(float(temp[0]),
                                                      n_subjects - 1)
    elif len(temp) == 2:
        thresh_clust = dict(start=float(temp[0]), step=float(temp[1]))
    else:
        raise ValueError()
    # To use an algorithm optimized for spatio-temporal clustering, we
    # just pass the spatial connectivity matrix (instead of spatio-temporal)
    print('Computing connectivity.')
    connectivity = mne.spatial_src_connectivity(src_fs)
    # Note that X needs to be a multi-dimensional array of shape
    # samples (subjects) x time x space, so we permute dimensions
    X = np.transpose(X, [2, 1, 0])
    print('Clustering.')
    stat_corr = 'cluster'
    if stat_corr == 'maxstat':
        # Permutation test with maximum statistic (conservative)
        T_obs, p_values, H0 = permutation_t_test(np.squeeze(X),
                                                 n_permutations=10000, tail=0,
                                                 n_jobs=N_JOBS, verbose=True)
    elif stat_corr == 'cluster':
        T_obs, clusters, cluster_p_values, H0 = clu = \
            spatio_temporal_cluster_1samp_test(
                X, connectivity=connectivity, n_jobs=N_JOBS,
                threshold=thresh_clust, buffer_size=None,
                n_permutations=5000, verbose=True, spatial_exclude=exclude,
                tail=tail)
    # Save observed T-values
    stc_T_obs = mne.SourceEstimate(np.transpose(T_obs),
                                   vertices=stc_cond1_fs.vertices, tmin=0,
                                   tstep=1, subject='fsaverage')
    if cond3:
        contrast_str = '{}-{}-{}{}_-_{}-{}-{}{}_vs_{}-{}-{}{}_-_{}-{}-{}{}'.format(
            cond1_mod, cond1_str, cond1_perm, prepend_if_not_empty(cond1_match_str, '-'),
            cond2_mod, cond2_str, cond2_perm, prepend_if_not_empty(cond2_match_str, '-'),
            cond3_mod, cond3_str, cond3_perm, prepend_if_not_empty(cond3_match_str, '-'),
            cond4_mod, cond4_str, cond4_perm, prepend_if_not_empty(cond4_match_str, '-'))
    else:
        contrast_str = '{}-{}-{}{}_vs_{}-{}-{}{}'.format(
            cond1_mod, cond1_str, cond1_perm, prepend_if_not_empty(cond1_match_str, '-'),
            cond2_mod, cond2_str, cond2_perm, prepend_if_not_empty(cond2_match_str, '-'))
    stc_fname = 'group_src-{}-tval_{}_{}_{}'.format(
        cond1_coh_str, contrast_str, freq_str, channel_type)
    stc_T_obs.save(op.join(dest_dir, stc_fname))

    if stat_corr == 'maxstat':
        stc_pval = mne.SourceEstimate(np.transpose(1-p_values),
                                      vertices=stc_cond1_fs.vertices, tmin=0,
                                      tstep=1, subject='fsaverage')
        stc_fname = 'group_src-{}-pval_{}_{}_{}'.format(
            cond1_coh_str, contrast_str, freq_str, channel_type)
        stc_pval.save(op.join(dest_dir, stc_fname))
    elif stat_corr == 'cluster':
        # Visualize the clusters
        good_cluster_inds = np.where(cluster_p_values < p_stat)[0]
        good_cluster_p = cluster_p_values[good_cluster_inds]

        # temp = np.zeros_like(np.transpose(T_obs))
        # for i, cl_idx in enumerate(good_cluster_inds):
        #     temp[clusters[cl_idx][1]] = 1-good_cluster_p[i]
        # stc_cluster = mne.SourceEstimate(temp, vertices=stc_cond1_fs.vertices,
        #                                 tmin=0, tstep=1, subject='fsaverage')
        # Now let's build a convenient representation of each cluster, where each
        # cluster becomes a "time point" in the SourceEstimate
        stc_cluster = summarize_clusters_stc(clu, vertices=fsave_vertices,
                                             subject='fsaverage',
                                             p_thresh=p_stat)
        if not mask:
            stc_fname = 'group_src-{}-clust_{}_{}_{}'.format(
                cond1_coh_str, contrast_str, freq_str, channel_type)
        else:
            id1 = mask.find('mask_')
            id2 = mask.find('.stc')
            mask_str = mask[id1+5:id2]
            stc_fname = 'group_src-{}-clust_{}_mask-{}_{}_{}'.format(
                cond1_coh_str, contrast_str, mask_str, freq_str, channel_type)

        stc_cluster.save(op.join(dest_dir, stc_fname))
        # Save cluster info as mat file
        clust_dict = {'T_obs': T_obs, 'clusters': clusters,
                      'cluster_p_values': cluster_p_values,
                      'H0': H0, 'good_cluster_inds': good_cluster_inds,
                      'good_cluster_p': good_cluster_p}
        savemat(op.join(dest_dir, '{}.mat'.format(stc_fname)), clust_dict)

elif any(stat == x for x in ['avg', 'max_vert']):
    X_cond1 = X_cond1.squeeze()
    X_cond2 = X_cond2.squeeze()
    if cond3:
        X_cond3 = X_cond3.squeeze()
        X_cond4 = X_cond4.squeeze()
    # Computing subject specific mean activities within the ROI. If not
    # specified use all vertices
    if mask and not anat_label:
        if stat == 'avg':
            m_cond1 = np.mean(X_cond1[include, :], 0).transpose()
            m_cond2 = np.mean(X_cond2[include, :], 0).transpose()
            if cond3:
                m_cond3 = np.mean(X_cond3[include, :], 0).transpose()
                m_cond4 = np.mean(X_cond4[include, :], 0).transpose()
        else:
            m_cond1 = np.amax(X_cond1[include, :], 0).transpose()
            m_cond2 = np.amax(X_cond2[include, :], 0).transpose()
            if cond3:
                m_cond3 = np.amax(X_cond3[include, :], 0).transpose()
                m_cond4 = np.amax(X_cond4[include, :], 0).transpose()
    else:
        if stat == 'avg':
            m_cond1 = np.mean(X_cond1, 0).transpose()
            m_cond2 = np.mean(X_cond2, 0).transpose()
            if cond3:
                m_cond3 = np.mean(X_cond3, 0).transpose()
                m_cond4 = np.mean(X_cond4, 0).transpose()
        else:
            m_cond1 = np.amax(X_cond1, 0).transpose()
            m_cond2 = np.amax(X_cond2, 0).transpose()
            if cond3:
                m_cond3 = np.amax(X_cond3, 0).transpose()
                m_cond4 = np.amax(X_cond4, 0).transpose()

    # Arranging data in a dataframe for plotting with seaborn
    sub = np.tile(np.array(range(1, n_subjects+1)), (n_conds))
    cond = np.tile('cond1', (n_subjects)).tolist()
    cond.extend(np.tile('cond2', (n_subjects)).tolist())
    stim_mod = np.tile(cond1_str, (n_subjects)).tolist()
    stim_mod.extend(np.tile(cond2_str, (n_subjects)).tolist())
    mod = np.tile(cond1_mod, (n_subjects)).tolist()
    mod.extend(np.tile(cond2_mod, (n_subjects)).tolist())
    perm = np.tile(cond1_perm, (n_subjects)).tolist()
    perm.extend(np.tile(cond2_perm, (n_subjects)).tolist())
    coh = np.concatenate((m_cond1, m_cond2), 0)
    if cond3:
        cond.extend(np.tile('cond3', (n_subjects)).tolist())
        cond.extend(np.tile('cond4', (n_subjects)).tolist())
        stim_mod.extend(np.tile(cond3_str, (n_subjects)).tolist())
        stim_mod.extend(np.tile(cond4_str, (n_subjects)).tolist())
        mod.extend(np.tile(cond3_mod, (n_subjects)).tolist())
        mod.extend(np.tile(cond4_mod, (n_subjects)).tolist())
        perm.extend(np.tile(cond3_perm, (n_subjects)).tolist())
        perm.extend(np.tile(cond4_perm, (n_subjects)).tolist())
        coh = np.concatenate((coh, m_cond3, m_cond4), 0)
    df = pd.DataFrame({'subID': pd.Categorical(sub),
                       'condition': pd.Categorical(cond),
                       'stim_modality': pd.Categorical(stim_mod),
                       'modality': pd.Categorical(mod),
                       'perm': pd.Categorical(perm),
                       'coherence': coh})

    # Do statistics
    # 1 sample ttest if only two conditions
    if not cond3:
        st, pval = stats.ttest_1samp(m_cond1-m_cond2, 0)
        # Correcting p-value if the test is one-sided
        if tail != 0:
            pval /= 2
    else:
        # 2x2 repeated measures ANOVA if 4 conditions
        aovrm = AnovaRM(df, 'coherence', 'subID', within=fact)
        res = aovrm.fit()
        print(res)

    # Plotting data
    if not cond3:
        pal = sns.color_palette()
        if cond1_mod == cond2_mod == 'aud':
            col = pal[0]
        elif cond1_mod == cond2_mod == 'vis':
            col = pal[1]
        else:
            col = None
        xticks = [cond1.split('-', maxsplit=1)[-1],
                  cond2.split('-', maxsplit=1)[-1]]
        if cond1_mod != cond2_mod and cond1_mod == 'vis':
            order = ['cond2', 'cond1']
            xticks.reverse()
        else:
            order = None
        # Plot data using seaborn
        sns.set_theme(style="ticks")
        f, ax = plt.subplots(figsize=(7, 6))
        # Plot the data with boxplot
        sns.boxplot(x='condition', y='coherence', data=df,
                    whis=[0, 100], color=col, width=.6, order=order)

        # Add in points to show each observation
        sns.stripplot(x='condition', y='coherence', data=df,
                      color=".3", linewidth=0, order=order)
        # Add t-test results as text
        plt.text(1, 1, 't-test:\nt({:d})={:.2f}\np({:d})={:.5f}'.format(
            n_subjects-1, st, tail, pval), horizontalalignment='right',
            verticalalignment='top', transform=ax.transAxes)
        # Tweak the visual presentation
        ax.set(xlabel='')
        sns.despine(trim=True, bottom=True)
        plt.xticks((0, 1), xticks)
        # Set title
        if mask:
            plt.title('{}, {}\n{}, hemi: {}'.format(
                cond1_coh_str, stat, mask, hemi))
        else:
            plt.title('{}, {}\nwhole brain'.format(cond1_coh_str, stat))
        plt.show()
    else:
        xticks = [cond1.split('-', maxsplit=2)[-1],
                  cond3.split('-', maxsplit=2)[-1]]
        order = None
        col = None
        # Plot data using seaborn
        sns.set_theme(style="ticks")
        f, ax = plt.subplots(figsize=(7, 5))
        # Plot the data with boxplot
        sns.boxplot(x=fact[0], y='coherence', hue=fact[1], data=df,
                    color=col, width=.6, order=order)

        # Add in points to show each observation
        # sns.stripplot(x=fact[0], y='coherence', hue=fact[1], data=df,
        #               color=".3", linewidth=0, order=order)
        # Add t-test results as text
        # plt.text(1, 1, 't-test:\nt({:d})={:.2f}\np({:d})={:.5f}'.format(
        #     n_subjects-1, st, tail, pval), horizontalalignment='right',
        #     verticalalignment='top', transform=ax.transAxes)
        # Tweak the visual presentation
        # ax.set(xlabel='')
        sns.despine(trim=True, bottom=True)
        # plt.xticks(tuple(range(4)), xticks)
        # Set title
        if mask:
            plt.title('{}, {}\n{}, hemi: {}'.format(
                cond1_coh_str, stat, mask, hemi))
        else:
            plt.title('{}, {}\nwhole brain'.format(cond1_coh_str, stat))
        plt.show()
