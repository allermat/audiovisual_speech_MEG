[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6822557.svg)](https://doi.org/10.5281/zenodo.6822557)

# Source code for research article: <br/>"Differential auditory and visual phase-locking are observed during audio-visual benefit and silent lip-reading for speech perception" 
**by Máté Aller, Heidi Solberg Økland, Lucy J. MacGregor, Helen Blank, and Matthew H. Davis**

Published in [Journal of Neuroscience](https://www.jneurosci.org/content/early/2022/06/24/JNEUROSCI.2476-21.2022)

Cite: Aller, M., Solberg Økland, H., MacGregor, L. J., Blank, H., & Davis, M. H. (2022). Differential auditory and visual phase-locking are observed during audio-visual benefit and silent lip-reading for speech perception. Journal of Neuroscience. https://doi.org/10.1523/JNEUROSCI.2476-21.2022



## Installation
1. Clone or download this repository
2. Download all dependencies into their individual subfolders within audiovisual_speech_MEG/code/toolbox/
3. Download the data folder from [here](https://osf.io/st6fe/) into the project root folder (audiovisual_speech_MEG/)
4. Build conda environment for python by executing:
	```
	conda env create -f environment.yml
	```

## Usage with MATLAB
1. Navigate to the project root folder (audiovisual_speech_MEG/) and edit AVSM_setupdir.m adding the path to the project root folder on your computer. Then execute AVSM_setupdir. This takes care of setting up the environment. 
2. Execute the matlab scripts found in results/ to reproduce behavioural and sensor space MEG results. 
3. All analysis code can be found in audiovisual_speech_MEG/code/analysis/

## Usage with Python
1. Setup environment by executing the commands below
	```
	cd path/to/project/folder
	conda activate avsm
	setenv MESA_GL_VERSION_OVERRIDE 3.3
	setenv PYTHONPATH `pwd`
	```
2. Update the project_dir variable in audiovisual_speech_MEG/code/analysis/megcoherence_utils.py to the project folder path on your system. 
3. To reproduce the source-space MEG results run the  audiovisual_speech_MEG/results/megcoherence_anatomical_roi_analysis.ipynb notebook in jupyter notebook. 
4. Source analysis pipeline can be found at audiovisual_speech_MEG/code/analysis/megcoherence_source_analysis_pipeline.txt. To reproduce source coherence maps run the corresponding line of code from the text file in an Ipython console. 

## Dependencies (MATLAB)
- [barwitherr](https://uk.mathworks.com/matlabcentral/fileexchange?q=barwitherr)
- [bayesFactor](https://github.com/klabhub/bayesFactor)
- [Colormaps](https://uk.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
- [crameri](https://uk.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps?s_tid=srchtitle)
- [fieldtrip](https://www.fieldtriptoolbox.org/download/)
- [gramm](https://github.com/piermorel/gramm)
- [randpermfull](https://uk.mathworks.com/matlabcentral/fileexchange/30189-randpermfull?s_tid=srchtitle)
- [shadedErrorBar](https://uk.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar)
- [suplabel](https://uk.mathworks.com/matlabcentral/fileexchange/7772-suplabel)

## Dependencies (Python)
All python dependencies are listed in audiovisual_speech_MEG/environment.yml

## Contributors
The MEG analysis code was written by Máté Aller. Preprocessing scripts (MEG, behavioural) were written by Heidi Solberg Økland. 
