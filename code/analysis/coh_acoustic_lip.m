function coh_acoustic_lip(varargin)
% Speech-speech coherence (lip aperture and acoustic envelope)
%
% This script uses FieldTrip to calculate coherence between the speech 
% envelope and lip aperture for all 275 experimental stimuli. We
% also compute coherence for mismatching lip aperture and speech envelopes
% to get a baseline for the coherence we could expect to see by chance.%
% 
% Based on Heidi Solberg Økland's script from February 2018
% 
% Copyright(C) Mate Aller 2020
% allermat@gmail.com

% Parsing input, checking matlab
p = inputParser;

validSentenceLength = {'whole','truncated'};

addOptional(p,'sentenceLength','truncated',@(x) ismember(x,validSentenceLength));
addOptional(p,'nPerm',5000,@(x) ismember(x,validSentenceLength));

parse(p,varargin{:});

sentenceLength= p.Results.sentenceLength;
nPerm = p.Results.nPerm;

% Load data
fprintf('Loading data...\n');
if strcmp(sentenceLength,'whole')
    % These are the whole sentences generated directly from the envelope
    % files scratch_get_envelopes_whole_sentences.m
    sourceDir = AVSM_setupdir('analysis_megcoherence_group');
    data = load(fullfile(sourceDir,'ft_aud_lip_env_whole.mat'));
    data = data.ftData;
else
    % The MEG datasets contain only truncated sentences
    sourceDir = AVSM_setupdir('analysis_megcoherence_sub','sub-01');
    data = load(fullfile(sourceDir,'trials_with_env_MEG'));
    data = data.trials_with_env_MEG;
end
destDir = AVSM_setupdir('analysis_megcoherence_group');

% --- 1) Frequency analysis ---

freqMin = 0.5; % Minimum frequency in Hz
freqMax = 20; % Maximum frequency in Hz
freqStep = 0.5; % Frequency steps in Hz
Padding = 10; % Total length of signal wanted in seconds
Smoothing = 0.5; % Spectral smoothing in Hz

fprintf('Doing Fourier analysis for frequencies between %d and %d...\n',freqMin,freqMax);

cfg = [];
cfg.output     = 'fourier'; % use Fourier analysis
cfg.method     = 'mtmfft'; % spectral smoothing using multitapers
cfg.foi        = freqMin:freqStep:freqMax; % from x to y Hz in steps of
cfg.pad        = Padding;
cfg.tapsmofrq  = Smoothing; % Hz - determines degree of smoothing
cfg.taper      = 'hanning';
cfg.keeptrials = 'yes';
cfg.channel    = {'AudEnv' 'LipEnv'};

fprintf('Doing Fourier analysis for frequencies between %d and %d...\n', freqMin, freqMax)
ftData_fourier    = ft_freqanalysis(cfg, data);

% --- 2) Coherence analysis, auditory envelope ---

cfg            = [];
cfg.method        = 'coh';
cfg.channelcmb    = {'LipEnv' 'AudEnv'};

% Do coherence analysis
fprintf('Computing coherence...\n');
coh = ft_connectivityanalysis(cfg, ftData_fourier);

cohPerm = cell(nPerm,1);
parfor iPerm = 1:nPerm
    ftData_temp = ftData_fourier;
    temp = ftData_temp.fourierspctrm;
    % Full permutation (a.k.a., derangement) can't be used here as this
    % data will be used to establish a null distribution which must be
    % sampled from all possible permutations. 
    temp(:,1,:) = temp(randperm(size(temp,1)),1,:);
    ftData_temp.fourierspctrm = temp;
    cohTemp = ft_connectivityanalysis(cfg, ftData_temp);
    % To spare some space I only save the relevant part of the cohspctrm 
    % field. The cohspctrm field contains a 2x2 matrix for each frequency
    % (i.e., 2x2xnFreq). From each frequency we just need one off-diagonal
    % element which represents the coherence between the two channels at
    % that frequency. The two off-diagonal entries are the same as
    % coherence is a symmetric measure. Each cell entry in cohPerm will 
    % contain a vector of length equal to the number of frequencies. 
    cohPerm{iPerm} = squeeze(cohTemp.cohspctrm(2,1,:));
end
% For saving, I concatenate the permuted coherence spectra into a single
% matrix of size nFrequencies x nPermutations. 
cohPerm = cat(2,cohPerm{:});
cohPerm = squeeze(cohPerm(2,1,:,:));
save(fullfile(destDir,sprintf('ftmeg_coh_acoustic_lip_%s',sentenceLength)),...
    'coh','cohPerm','-v7.3');

% Computing the powerspectra of the stimuli
cfg = [];
cfg.output     = 'pow'; % use Fourier analysis
cfg.method     = 'mtmfft'; % spectral smoothing using multitapers
cfg.foi        = freqMin:freqStep:freqMax; % from x to y Hz in steps of
cfg.pad        = Padding;
cfg.tapsmofrq  = Smoothing; % Hz - determines degree of smoothing
cfg.taper      = 'hanning';
cfg.keeptrials = 'yes';
cfg.channel    = {'AudEnv' 'LipEnv'};
ftData_pow    = ft_freqanalysis(cfg, data);

% Computing the powerspectra of the baseline
cfg = [];
cfg.latency = [0,1]; % Based on visual inspection this period baseline in all trials
ftData_baseline = ft_selectdata(cfg,data);

cfg = [];
cfg.output     = 'pow'; % use Fourier analysis
cfg.method     = 'mtmfft'; % spectral smoothing using multitapers
cfg.foi        = freqMin:freqStep:freqMax; % from x to y Hz in steps of
cfg.pad        = Padding;
cfg.tapsmofrq  = Smoothing; % Hz - determines degree of smoothing
cfg.taper      = 'hanning';
cfg.keeptrials = 'yes';
cfg.channel    = {'AudEnv' 'LipEnv'};
ftData_pow_baseline    = ft_freqanalysis(cfg, ftData_baseline);

save(fullfile(destDir,sprintf('ftmeg_coh_acoustic_lip_%s',sentenceLength)),...
    'ftData_pow','ftData_pow_baseline','-append');
end
