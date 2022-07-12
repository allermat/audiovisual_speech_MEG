# -*- coding: utf-8 -*-
"""Plot source space data

    Parameters
    ----------
    file_name : str
        Name of the file to be plotted
    clims : str
        Color bar limits as a hyphen separated list
    hemi : str
        Hemisphere to plot (both, split, lh or rh)
        Default: both
    surface : str
        Type of surface (inflated, inflated_pre, pial_semi_inflated)
        Default: inflated
    source_space : str
        Source space type (surf of vol) Not used currently
        Default: surf
    time_viewer : str
        Time viewer on or off
    views : str
        Views for plotting  [lat, med, ros, cau, dor, ven, fro, par]
        Default: lat

Created on Fri May 22 11:03:35 2020

@author: Máté Aller
"""

import os.path as op
from os import listdir
import mne
import numpy as np
from argparse import ArgumentParser
from scipy import stats as stats
from cmcrameri import cm
from megcoherence_utils import project_dir
# import re
# from matplotlib import pyplot as plt
# from matplotlib import cm
# from matplotlib.colors import ListedColormap

# Parse input
parser = ArgumentParser()
parser.add_argument('file_name', help='Name of the file to be plotted')
parser.add_argument('-cl', '--clims', help='Color bar limits',
                    default='')
parser.add_argument('-bg', '--background', help='Background colour of source plot (matplotlib color)',
                    default='k')
parser.add_argument('-ct', '--channel_type', help='Channel type (mag or grad)',
                    default='mag')
parser.add_argument('-cmap', '--color_map', help='Color map to use (auto or any from cmcrameri)',
                    default='auto')
parser.add_argument('-he', '--hemi', help='Hemisphere to plot (both, split, lh or rh)',
                    default='both')
parser.add_argument('-sd', '--subdir', help='Subdirectory to dave data',
                    default='')
parser.add_argument('-su', '--surface', help='Type of surface (inflated, inflated_pre, pial_semi_inflated)',
                    default='inflated_pre')
parser.add_argument('-ss', '--source_space', help='Source space type (surf or vol)',
                    default='surf')
parser.add_argument('-tv', '--time_viewer', help='Time viewer on or off',
                    default='on')
parser.add_argument('-v', '--views', help='Views for plotting  [lat, med, ros, cau, dor, ven, fro, par]',
                    default='lat')
parser.add_argument('-l', '--label', help='Anatomical label to display',
                    default='')
parser.add_argument('-fl', '--func_label_fname', help='Functional label file',
                    default='')
parser.add_argument('-flc', '--func_label_color', help='Functional label color',
                    default='b')
parser.add_argument('-flm', '--func_label_mode', help='Functional label display mode',
                    default='draw')

args = parser.parse_args()
file_name = args.file_name
background = args.background
clims = args.clims
channel_type = args.channel_type
color_map = args.color_map
hemi = args.hemi
subdirectory = args.subdir
surface = args.surface
source_space = args.source_space
time_viewer = args.time_viewer
views = args.views.split('-')
label = args.label
func_label_fname = args.func_label_fname
func_label_color = args.func_label_color
func_label_mode = args.func_label_mode

mne.viz.set_3d_backend('pyvista')

# Find subID
subID = file_name.split('_')[0]
is_corr = file_name.split('_')[1] == 'src-corr'

if time_viewer == 'on':
    time_viewer = True
elif time_viewer == 'off':
    time_viewer = False
else:
    raise ValueError()

assert any(func_label_mode == x for x in ['draw', 'mask'])
if func_label_mode == 'mask':
    assert func_label_fname # should not be empty

subjects_dir = op.join(project_dir, 'data', 'derivatives', 'anat')
if subdirectory:
    dest_dir = op.join(project_dir, 'data', 'derivatives', 'megcoherence', subID, subdirectory)
else:
    dest_dir = op.join(project_dir, 'data', 'derivatives', 'megcoherence', subID)
subj_dir_list = [f for f in listdir(subjects_dir)
                 if op.isdir(op.join(subjects_dir, f))]
n_subj = sum(1 for s in subj_dir_list if s.startswith('sub'))

if label:
    if label == 'Occ_cortex':
        label = ['G_and_S_occipital_inf', 'G_cuneus', 'G_occipital_middle',
                 'G_occipital_sup', 'Pole_occipital',
                 'S_oc_middle_and_Lunatus', 'S_oc_sup_and_transversal']
    elif label == 'STG':
        label = ['G_temp_sup-Lateral']
    else:
        label = [label]

    # Load label
    def fun_read_label(lab, h):
        return mne.read_labels_from_annot('fsaverage', parc='aparc.a2009s',
                                          hemi=h, subjects_dir=subjects_dir,
                                          regexp=lab)[0]
    if any(hemi == s for s in ['both', 'split']):
        temp = [fun_read_label(s, 'lh') for s in label]
        anat_label_lh = temp[0]
        for x in temp[1:]:
            anat_label_lh += x
        temp = [fun_read_label(s, 'rh') for s in label]
        anat_label_rh = temp[0]
        for x in temp[1:]:
            anat_label_rh += x
    else:
        temp = [fun_read_label(s, hemi) for s in label]
        anat_label = temp[0]
        for x in temp[1:]:
            anat_label += x

# Load mask file for displaying outline
if func_label_fname:
    func_label_stc = mne.read_source_estimate(op.join(dest_dir, func_label_fname))
    # func_label_stc = func_label_stc.crop(0,0)
    fs_src_fname = subjects_dir + '/fsaverage/bem/fsaverage-ico-5-src.fif'
    src_fs = mne.read_source_spaces(fs_src_fname)
    func_label = mne.stc_to_label(func_label_stc, src=src_fs, smooth=True,
                                  connected=False, subjects_dir=subjects_dir)

# Set colorbar limits
if not clims:
    clims = 'auto'
else:
    temp = clims.split('-')
    if len(temp) == 3:
        clims = dict(kind='percent', pos_lims=list(map(float, temp)))
    elif len(temp) == 5:
        kind = temp[0]
        lim = temp[1]
        if kind == 'pval':
            assert subID == 'group'
            # Compute t values corresponding to the given pvalues
            if lim == 'pos_lims':
                lims = [-stats.distributions.t.ppf(float(t)/2., n_subj-1)
                        for t in temp[2:]]
            elif lim == 'lims':
                lims = [-stats.distributions.t.ppf(float(t), n_subj-1)
                        for t in temp[2:]]
            # Converting critical t value to R if the data are source
            # correlations
            if is_corr:
                lims = [t / np.sqrt(n_subj - 2 + t**2) for t in lims]
        elif kind == 'value':
            lims = temp[2:]

        if lim == 'pos_lims':
            clims = dict(kind='value', pos_lims=list(map(float, lims)))
        elif lim == 'lims':
            clims = dict(kind='value', lims=list(map(float, lims)))
        else:
            raise ValueError()
    else:
        raise ValueError()

if color_map == 'auto':
    cmap = 'auto'
else:
    cmap = getattr(cm, color_map)

# Load source estimates
stc = mne.read_source_estimate(op.join(dest_dir, file_name))

if func_label_mode == 'mask':
    temp_data = stc.data
    temp_data[func_label_stc.data[:, 0] == 0] = 0
    stc.data = temp_data
    # stc = stc.in_label(func_label[0])

# Set color map based on file name - not currently used
# from IPython.core.debugger import Pdb; ipdb = Pdb(); ipdb.set_trace()
# p = re.compile('[a-zA-Z0-9-]+_([a-zA-Z0-9-]+)_.*')
# m = p.match(stc_fname)
# if m.group(1) == 'src-coh-tval':
#     cmap = 'mne'
# elif m.group(1) == 'src-coh-clust':
#     cmap = 'auto'
#     # hot = cm.get_cmap('hot',256)
#     # cmap = ListedColormap(np.flipud(hot(range(256))))
# else:
#     cmap = 'auto'

# Visualize the reconstructed source activity
if subID == 'group':
    brain = stc.plot(hemi=hemi, subjects_dir=subjects_dir, surface=surface,
                     subject='fsaverage', views=views,
                     clim=clims, colormap=cmap, time_viewer=time_viewer,
                     background=background)
else:
    brain = stc.plot(hemi=hemi, subjects_dir=subjects_dir, surface=surface,
                     subject=subID, views=views,
                     clim=clims, colormap=cmap, time_viewer=time_viewer,
                     background=background)
if label:
    if any(hemi == s for s in ['both', 'split']):
        brain.add_label(anat_label_lh, borders=True, color='k')
        brain.add_label(anat_label_rh, borders=True, color='k')
    else:
        brain.add_label(anat_label, borders=True, color='k')

if func_label_fname and func_label_mode == 'draw':
    # if any(hemi == s for s in ['both', 'split']):
    #     brain.add_label(anat_label_lh, borders=True, color='k')
    #     brain.add_label(anat_label_rh, borders=True, color='k')
    # else:
    if func_label[0] is not None:
        brain.add_label(func_label[0], borders=True, color=func_label_color)
    if func_label[1] is not None:
        brain.add_label(func_label[1], borders=True, color=func_label_color)
        
