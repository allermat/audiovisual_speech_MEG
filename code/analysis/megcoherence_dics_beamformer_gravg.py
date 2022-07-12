# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:03:35 2020

@author: ma09
"""
import os.path as op
import mne
from argparse import ArgumentParser
from megcoherence_utils import project_dir

parser = ArgumentParser()
parser.add_argument('-ct', '--channel_type', help='Channel type (mag or grad)',
                    default='mag')
parser.add_argument('-m', '--modality', help='Sensory modality (A or V)',
                    default='A')
parser.add_argument('-o', '--output', help='Output type (true or perm)',
                    default='true')
parser.add_argument('-ss', '--source_space', help='Source space type (surf or vol)',
                    default='surf')
args = parser.parse_args()
channel_type = args.channel_type
modality = args.modality
output = args.output
source_space = args.source_space

if modality == 'A':
    cond_str = 'allAud'
    mod_str = 'aud'
else:
    cond_str = 'allVis'
    mod_str = 'vis'
freq_str = '{}-{}Hz'.format(2,6)

n_subjects = 14
subIDlist = ['sub-' + str(n).zfill(2) for n in list(range(1,n_subjects+1))]

subjects_dir = op.join(project_dir, 'data', 'derivatives', 'anat')
dest_dir = op.join(project_dir, 'data', 'derivatives', 'megcoherence', 'group') 
# load fsaverage source space for morphing
mne.datasets.fetch_fsaverage(subjects_dir)  # ensure fsaverage src exists
fname_fs_src = subjects_dir + '/fsaverage/bem/fsaverage-vol-5-src.fif'
src_fs = mne.read_source_spaces(fname_fs_src)

stc_morphed = []
for subID in subIDlist:
    
    source_dir = op.join(project_dir, 'data', 'derivatives', 'megcoherence', subID)    
    # Load precomputed forward solution and source estimates
    # ------------------------
    if source_space == 'vol':
        forward_fname = op.join(source_dir,''.join([subID,'-fwd.fif']))
        source_fname = op.join(source_dir,''.join([subID,'_',modality,'_stc-vl.stc']))
        
        forward = mne.read_forward_solution(forward_fname)
        stc = mne.read_source_estimate(source_fname)
        
        morph = mne.compute_source_morph(
            forward['src'], subject_from=subID, src_to=src_fs,
            subjects_dir=subjects_dir,
            niter_sdr=[10, 10, 5], niter_affine=[10, 10, 5],  # just for speed
            verbose=True)
        
    else:
        forward_fname = op.join(source_dir,''.join([subID,'_surf-fwd.fif']))
        forward = mne.read_forward_solution(forward_fname)

        stc_fname = '{}_mne_src-coh_{}-{}_{}_{}_{}-lh.stc'.format(
            subID,mod_str,cond_str,freq_str,channel_type,output)
        stc = mne.read_source_estimate(op.join(source_dir,stc_fname))
        # Morph source estimates to fsaverage for group analysis
        morph = mne.compute_source_morph(
            forward['src'], subject_from=subID, subject_to='fsaverage',
            subjects_dir=subjects_dir)
    
    stc_fs = morph.apply(stc)
    # stc_fs.save(op.join(source_dir,''.join([subID,'_',modality,'_fsavg_stc'])))
    stc_morphed.append(stc_fs)
    
# Average source estimates across subjects
stc_gravg = stc_fs.copy()
for subIdx in range(0,n_subjects):
    stc_gravg._data += stc_morphed[subIdx].data
stc_gravg._data /= n_subjects

# Save group average source estimate
if source_space == 'vol':
    gravg_fname = op.join(dest_dir,''.join(['group_',modality,'_stc']))
else:
    gravg_fname = stc_fname = '{}_mne_src-coh_{}-{}_{}_{}_{}'.format(
            'group',mod_str,cond_str,freq_str,channel_type,output)
stc_gravg.save(op.join(dest_dir,gravg_fname))