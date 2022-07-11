# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:03:35 2020

@author: ma09
"""
from argparse import ArgumentParser
import os.path as op
import re
import mne
import numpy as np
from scipy import stats as stats

parser = ArgumentParser()

parser.add_argument('timg', help='T-image file name')
parser.add_argument('-m2', '--mask2', help='Second T-image file name',
                    default='')
parser.add_argument('-ps', '--p_stat', help='Statistical comparison p threshold',
                    default='0.05')
parser.add_argument('-t', '--tail', help='Statistical comparison tail (-1, 0, 1)',
                    default='1')
parser.add_argument('-mfn', '--mask_fname', help='Mask file name for saving',
                    default='default')
args = parser.parse_args()
timg = args.timg
mask2 = args.mask2
p_stat = float(args.p_stat)
tail = int(args.tail)
mask_fname = args.mask_fname

n_subjects = 14


subjects_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects',
                       'AVSpeechMEG', 'data', 'derivatives', 'anat')
dest_dir = op.join('/imaging', 'davis', 'users', 'ma09', 'Projects',
                   'AVSpeechMEG', 'data', 'derivatives', 'megcoherence',
                   'group')
# Load mask image(s)
stc_timg = mne.read_source_estimate(op.join(dest_dir,timg))
if mask2:
    stc_mask2 = mne.read_source_estimate(op.join(dest_dir,mask2))
    data_mask2 = stc_mask2.data

# check if the input is a t image or results of cluster statistics
m = re.search('-clust',timg)
if m:
    clust = stc_timg.data
    # the first column of cluster stats is all clusters together
    good_cluster = np.where(np.absolute(clust[:,0]) > 0)[0]
    # Creating mask data
    data = np.zeros_like(clust)
else:
    # Find above threshold T-values
    tval_timg = stc_timg.data
    if tail == 0:
        t_thresh = -stats.distributions.t.ppf(p_stat / 2., n_subjects - 1)
        good_cluster = tval_timg > t_thresh
    elif tail == 1:
        t_thresh = -stats.distributions.t.ppf(p_stat, n_subjects - 1)
        good_cluster = tval_timg > t_thresh
    elif tail == -1:
        t_thresh = -stats.distributions.t.ppf(p_stat, n_subjects - 1)
        good_cluster = tval_timg < -t_thresh
    # Creating mask data
    data = np.zeros_like(tval_timg)
# Applying mask
if not mask2:
    data[good_cluster] = 1
else:
    data[np.logical_and(good_cluster,data_mask2 != 0)] = 1

# Write mask as stc
stc_mask = mne.SourceEstimate(data, vertices=stc_timg.vertices, 
                              tmin=0, tstep=1, subject='fsaverage')
if m:
    stc_fname = 'group_src-mask_{}_{}'.format(mask_fname,'clust')
else:
    stc_fname = 'group_src-mask_{}_{}'.format(mask_fname,'tval_p{}'.format(str(p_stat)[2:]))
stc_mask.save(op.join(dest_dir,stc_fname))



