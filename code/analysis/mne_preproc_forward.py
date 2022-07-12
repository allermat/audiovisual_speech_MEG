"""
.. _tut-forward:

Head model and forward computation
==================================

The aim of this tutorial is to be a getting started for forward
computation.

For more extensive details and presentation of the general
concepts for forward modeling, see :ref:`ch_forward`.
"""

import os.path as op
import mne
from argparse import ArgumentParser
from megcoherence_utils import project_dir

parser = ArgumentParser()
parser.add_argument('fwd_type', '--ft', help='Type of forward model (vol or surf)',
                    default='surf')
parser.add_argument('source_space', '--ss', help='Type of forward model (indiv or fsaverage)',
                    default='indiv')

args = parser.parse_args()
fwd_type = args.fwd_type
source_space = args.source_space

# The paths to Freesurfer reconstructions
subjects_dir = op.join(project_dir, 'data', 'derivatives', 'anat')

subIDlist = ['sub-' + str(n).zfill(2) for n in list(range(1, 15))]

if source_space == 'fsaverage':
    # In this case the source space is defined on the fsaverage and morphed
    # into each individual subject's anatomy. 
    fsavg_src_fname = op.join(subjects_dir, 'fsaverage', 'fsaverage-src.fif')
    if op.isfile(fsavg_src_fname):
        fsaverage_src = mne.read_source_spaces(fsavg_src_fname)
    else:
        fsaverage_src = mne.setup_source_space('fsaverage', spacing='oct6',
                                               subjects_dir=subjects_dir,
                                               n_jobs=2, add_dist='patch')
        mne.write_source_spaces(fsavg_src_fname, fsaverage_src, overwrite=True)

for subID in subIDlist:
    raw_data_dir = op.join(project_dir, 'data', 'derivatives', 'maxfilter', subID)
    # the raw file containing the channel location + types
    raw_fname = op.join(raw_data_dir, 'block1_raw.fif')
    destDir = op.join(project_dir, 'data', 'derivatives', 'megcoherence', subID)
    # Visualizing the coregistration
    # ------------------------------
    #
    # The coregistration is the operation that allows to position the head and the
    # sensors in a common coordinate system. In the MNE software the transformation
    # to align the head and the sensors in stored in a so-called **trans file**.
    # It is a FIF file that ends with ``-trans.fif``. It can be obtained with
    # :func:`mne.gui.coregistration` (or its convenient command line
    # equivalent :ref:`gen_mne_coreg`), or mrilab if you're using a Neuromag
    # system.
    #
    # Here we assume the coregistration is done, so we just visually check the
    # alignment with the following code.
    
    # mne.gui.coregistration(subject=subID,inst=raw_fname)
    
    # The transformation file obtained by coregistration
    trans = op.join(raw_data_dir,''.join([subID,'-trans.fif']))
    
    info = mne.io.read_info(raw_fname)
    # Here we look at the dense head, which isn't used for BEM computations but
    # is useful for coregistration.
    # mne.viz.plot_alignment(info, trans, subject=subID, dig=True,
    #                        meg=['helmet', 'sensors'], subjects_dir=subjects_dir,
    #                        surfaces='head')
    
    # Compute Source Space
    # --------------------
    #
    # The source space defines the position and orientation of the candidate source
    # locations. There are two types of source spaces:
    #
    # - **surface-based** source space when the candidates are confined to a
    #   surface.
    #
    # - **volumetric or discrete** source space when the candidates are discrete,
    #   arbitrarily located source points bounded by the surface.
    #
    # **Surface-based** source space is computed using
    # :func:`mne.setup_source_space`, while **volumetric** source space is computed
    # using :func:`mne.setup_volume_source_space`.
    #
    
    # To compute a volume based source space defined with a grid of candidate
    # dipoles inside the brain (requires the :term:`BEM` surfaces) you can use the
    # following.
    
    if fwd_type == 'vol':
        surface = op.join(subjects_dir, subID, 'bem', 'inner_skull.surf')
        vol_src = mne.setup_volume_source_space(subID, subjects_dir=subjects_dir,
                                                surface=surface)
        print(vol_src)
    elif fwd_type == 'surf':
        if source_space == 'fsaverage':
            # Morph the source space to the current subject
            src = mne.morph_source_spaces(fsaverage_src, subID,
                                          subjects_dir=subjects_dir)
        else:
            src = mne.setup_source_space(subID, spacing='oct6',
                                         add_dist='patch',
                                         subjects_dir=subjects_dir)
            print(src)
    else:
        raise ValueError

    # mne.viz.plot_bem(subject=subID, subjects_dir=subjects_dir,
                     # brain_surfaces='white', src=vol_src, orientation='coronal')

    # Compute forward solution
    # ------------------------
    #
    # We can now compute the forward solution.
    # To reduce computation we'll just compute a single layer BEM (just inner
    # skull) that can then be used for MEG (not EEG).
    # We specify if we want a one-layer or a three-layer BEM using the
    # ``conductivity`` parameter.
    # The BEM solution requires a BEM model which describes the geometry
    # of the head the conductivities of the different tissues.

    conductivity = (0.3,)  # for single layer
    # conductivity = (0.3, 0.006, 0.3)  # for three layers
    model = mne.make_bem_model(subject=subID, ico=4,
                               conductivity=conductivity,
                               subjects_dir=subjects_dir)
    bem = mne.make_bem_solution(model)

    # Note that the :term:`BEM` does not involve any use of the trans file. The BEM
    # only depends on the head geometry and conductivities.
    # It is therefore independent from the MEG data and the head position.

    # Let's now compute the forward operator, commonly referred to as the
    # gain or leadfield matrix.
    # See :func:`mne.make_forward_solution` for details on the meaning of each
    # parameter.

    if fwd_type == 'vol':
        fwd = mne.make_forward_solution(raw_fname, trans=trans, src=vol_src, bem=bem,
                                        meg=True, eeg=False, mindist=5.0, n_jobs=2)
        fwd_fname = op.join(destDir, ''.join([subID, '_vol-fwd.fif']))
    elif fwd_type == 'surf':
        fwd = mne.make_forward_solution(raw_fname, trans=trans, src=src, bem=bem,
                                        meg=True, eeg=False, mindist=5.0, n_jobs=2)
        fwd_fname = op.join(destDir, ''.join([subID, '_surf-fwd.fif']))
    print(fwd)
    # Write forward solution to disk
    mne.write_forward_solution(fwd_fname, fwd, overwrite=True)
