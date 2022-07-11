function coh_sensor_aud_vis(subID,varargin)
% Coherence analysis for basic auditory/visual effects
%
%
% This script uses FieldTrip to calculate coherence between the speech 
% envelope + lip aperture and MEG. It collapses over all conditions where
% there is auditory speech (AO20, AO70, AV20 and AV70) and visual speech 
% (AV20, AV70 and VO). To ensure that no one condition has more influence 
% than another, we pick an equal number of trials from all relevant 
% conditions for a given subject before collapsing over auditory/visual 
% trials. We also make a permuted dataset with random pairings of speech 
% envelope + lip aperture with the MEG data to look at chance-level 
% coherence. The coherence is only calculated for gradiometers. 
% 
% Based on Heidi Solberg Økland's script from November2017-March 2018
% 
% Copyright(C) Mate Aller 2020
% allermat@gmail.com

% Parsing input, checking matlab
p = inputParser;

validMegChannels = {'MEG','MEGGRAD','MEGMAG'};
validCohAnalyses = {'all','basicAudVis'};

addRequired(p,'subID',@ischar);
addParameter(p,'nPerm',100,@(x) validateattributes(x,{'double'},{'scalar',...
                                'integer','positive'}));
addParameter(p,'chanMeg','MEG',@(x) ismember(x,validMegChannels));
addParameter(p,'combinePlanar',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'doDerangement',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'doSquaredCoh',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'cohAnalyses','all',@(x) ismember(x,validCohAnalyses));

parse(p,subID,varargin{:});

subID = p.Results.subID;
nPerm = p.Results.nPerm;
chanMeg = p.Results.chanMeg;
combinePlanar = p.Results.combinePlanar;
doDerangement = p.Results.doDerangement;
doSquaredCoh = p.Results.doSquaredCoh;
cohAnalyses = p.Results.cohAnalyses;

sourceDir = AVSM_setupdir('analysis_megcoherence_sub',subID);

% key: 16 = AO 20%, 32 = AO 70%, 48 = AV 20%, 64 = AV 70%, 80 = VO
allCondLab = {'AO20','AO70','AV20','AV70','VO'};
allTrigCodes = [16,32,48,64,80];

audCondLab = {'AO20','AO70','AV20','AV70'};
visCondLab = {'AV20','AV70','VO'};
audOnlyCondLab = {'AO20','AO70'};

fprintf('====== Doing auditory and visual coherence analysis for %s ======\n\n',subID)

% Load data
fprintf('Loading data...\n');
ftData = load(fullfile(sourceDir,'trials_with_env_MEG.mat'));
ftData = ftData.trials_with_env_MEG;

% 1) Find out how many good trials there are per condition
nTrialsCond = arrayfun(@(x) sum(ismember(ftData.trialinfo,x)),allTrigCodes);

% Find out the maximum number of trials to include per condition by
% finding the minimum number of good trials in one of the conditions
min_nTrials = min(nTrialsCond);
fprintf('\n-----> Number of trials per condition going into the analysis is %d. \n\n', min_nTrials)

% Find trial indices for each condition
trlIdxPerCond = cell(size(allCondLab));
for j = 1:length(allCondLab)
    % Find the trials that have a certain trigger code (condition)
    trlIdx = find(ftData.trialinfo==allTrigCodes(j));
    % Choose an equal amount of trials for all conditions (omitting the
    % highest numbers as subjects are usually tired towards the end..)
    trlIdxPerCond{j} = trlIdx(1:min_nTrials);
end

% Collect all auditory condition trial indices
trlIdxAllAud = cat(1,trlIdxPerCond{ismember(allCondLab,audCondLab)});
trlIdxAllAud = sort(trlIdxAllAud);

% Collect all visual condition trial indices
trlIdxAllVis = cat(1,trlIdxPerCond{ismember(allCondLab,visCondLab)});
trlIdxAllVis = sort(trlIdxAllVis);

% Collect all visual condition trial indices
trlIdxAOall = cat(1,trlIdxPerCond{ismember(allCondLab,audOnlyCondLab)});
trlIdxAOall = sort(trlIdxAOall);

% Combine gradiometers if necesarry
if combinePlanar
    ftData = ft_combineplanar([],ftData);
end

% FFT settings
freqMin = 1; % Minimum frequency in Hz
freqMax = 20; % Maximum frequency in Hz
freqStep = 1; % Frequency steps in Hz
Padding = 10; % Total length of signal wanted in seconds
Smoothing = 0.6; % Spectral smoothing in Hz

cfg = [];
cfg.output = 'fourier'; % use Fourier analysis
cfg.method = 'mtmfft'; % spectral smoothing using multitapers
cfg.foi = freqMin:freqStep:freqMax; % from x to y Hz in steps of
cfg.pad = Padding;
cfg.tapsmofrq  = Smoothing; % Hz - determines degree of smoothing
cfg.taper = 'hanning';
cfg.keeptrials = 'yes';
cfg.trials = 'all';

% Do Fourier analysis
fprintf('Doing Fourier analysis for frequencies between %d and %d...\n',freqMin,freqMax)
ftData_fourier = ft_freqanalysis(cfg,ftData);

% Check for NaNs
if any(isnan(ftData_fourier.fourierspctrm))
    error('The fourier transfor has NaNs in it. Something is wrong :(')
end

% Generate all envelope and condition combinations for the coherence
% analysis
if strcmp(cohAnalyses,'all')
    analysisConds = transpose(cat(2,allCondLab,{'AOall','allAud','allVis'}));
    analysisTrlIdx = transpose(cat(2,trlIdxPerCond,{trlIdxAOall,trlIdxAllAud,trlIdxAllVis}));
    analysisDetails = repmat(analysisConds,2,1);
    analysisDetails(:,2) = repmat(analysisTrlIdx,2,1);
    temp = repmat({'aud','vis'},numel(analysisConds),1);
    analysisDetails(:,3) = temp(:);
    % Remove combinations which are not needed (auditory coherence in allVis
    % and visual coherence in allAud)
    analysisDetails(ismember(analysisDetails(:,1),'allAud') & ismember(analysisDetails(:,3),'vis'),:) = [];
    analysisDetails(ismember(analysisDetails(:,1),'allVis') & ismember(analysisDetails(:,3),'aud'),:) = [];
elseif strcmp(cohAnalyses,'basicAudVis')
    analysisDetails = {'allAud',trlIdxAllAud,'aud';...
                       'allVis',trlIdxAllVis,'vis'};
else
    error('Invalid coherence analysis');
end

for i = 1:size(analysisDetails,1)
    
    if strcmp(analysisDetails{i,3},'aud')
        chanEnv = 'AudEnv';
    else
        chanEnv = 'LipEnv';
    end
    actTrlIdx = analysisDetails{i,2};
    % Compute true coherence
    cfg = [];
    cfg.method = 'coh';
    % This makes sure that the coherency is kept (i.e both average
    % direction and magnitude)
    cfg.complex = 'abs';
    cfg.channelcmb = {chanMeg,chanEnv};
    cfg.trials = actTrlIdx;
    fprintf('Computing coherence: %s-%s...\n',analysisDetails{i,3},analysisDetails{i,1})
    ftData_coh = ft_connectivityanalysis(cfg,ftData_fourier);
    % Fieldtrip computes the absolute value of coherency, see in 
    % ft_connectivity_corr.m and this thread 
    % (https://mailman.science.ru.nl/pipermail/fieldtrip/2014-February/020455.html) 
    % on the fieldtrip mailing list. If we want the magnitude-squared
    % coherence, we need to square the values. 
    if doSquaredCoh
        ftData_coh.cohspctrm = (ftData_coh.cohspctrm).^2;
    end
    
    % Save
    fprintf('Saving...\n')
    % Define output filenames
    fileName = [subID,'_ftmeg_coh_',analysisDetails{i,3},'-',analysisDetails{i,1}];
    if combinePlanar, fileName = strcat(fileName,'_combPlanar'); end
    save(fullfile(sourceDir,fileName),'ftData_coh','-v7.3')
    % Free up space
    ftData_coh = [];
    
    if ismember(analysisDetails{i,1},{'VO','AOall','allAud','allVis'})
        % Make permuted data (swap envelopes around) and compute average 
        % permuted coherence
        
        % Make sure that the random number generator starts from a different
        % seed for each subject/each time the script is run
        rng('shuffle');
        
        ftDataCell_coh_perm = cell(nPerm,1);
        % Loop over a pre-defined number of permutations to get average
        % permuted coherence
        fprintf('\nRunning permutations\n\n');
        for iPerm = 1:nPerm
            % Select the required trials
            cfg = [];
            cfg.trials = actTrlIdx;
            ftData_temp = ft_selectdata(cfg,ftData_fourier);
            % Shuffle envelopes
            temp = ftData_temp.fourierspctrm;
            envIdx = ismember(ftData_temp.label,chanEnv);
            % I apply exact permutations (i.e., derangements) which makes 
            % sure that no item falls in its respective original place 
            % after permutation. 
            if doDerangement
                temp(:,envIdx,:) = temp(randpermfull(size(temp,1)),envIdx,:);
            else
                temp(:,envIdx,:) = temp(randperm(size(temp,1)),envIdx,:);
            end
            ftData_temp.fourierspctrm = temp;
            % Compute permuted coherence
            cfg = [];
            cfg.method = 'coh';
            cfg.complex = 'abs';
            cfg.channelcmb = {chanMeg,chanEnv};
            temp_coh = ft_connectivityanalysis(cfg,ftData_temp);
            % If we want the magnitude-squared coherence, we need to square
            % the values here. See above for details at true coherence 
            % computation. 
            if doSquaredCoh
                temp_coh.cohspctrm = (temp_coh.cohspctrm).^2;
            end
            ftDataCell_coh_perm{iPerm} = temp_coh;
        end
        
        % I use ftData_coh as name as later I want to save the data with
        % this variable name. 
        ftData_coh = ftDataCell_coh_perm{1};
        temp = [ftDataCell_coh_perm{:}];
        ftData_coh.cohspctrm = mean(cat(4,temp.cohspctrm),4);
        
        % Save permuted coherence
        fprintf('Saving permuted coherence...\n');
        fileName = [subID,'_ftmeg_coh_',analysisDetails{i,3},'-',analysisDetails{i,1},'_perm'];
        if combinePlanar, fileName = strcat(fileName,'_combPlanar'); end
        save(fullfile(sourceDir,fileName),'ftData_coh','-v7.3');
        
    end
end

end


