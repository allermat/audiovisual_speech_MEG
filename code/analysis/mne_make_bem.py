#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 16:30:53 2020

@author: ma09
"""
import os.path as op
import mne
from megcoherence_utils import project_dir

subIDlist = ['sub-' + str(n).zfill(2) for n in list(range(1,15))]

# The paths to Freesurfer reconstructions
subjects_dir = op.join(project_dir, 'data', 'derivatives', 'anat')
for subID in subIDlist:
    mne.bem.make_watershed_bem(subID,subjects_dir=subjects_dir)
