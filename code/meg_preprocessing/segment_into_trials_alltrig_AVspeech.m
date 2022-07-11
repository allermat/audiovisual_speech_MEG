%% Segment data into trials
%
% Heidi Solberg Okland
%
% Use FieldTrip functions ft_definetrial and ft_preprocessing to segment
% the continuous recording into trials. Do this after ICA/adjusting head
% position.
% The script takes a few minutes to run on one dataset.

clear all; close all;

addpath(genpath('/imaging/local/software/fieldtrip/fieldtrip-latest'));

dataPath = '';

subjPaths = {};
    
fileName_in = 'concatenated_icaed_raw_trans.fif';

triggerCode = [16 32 48 64 80];
% key: 16 = AO 20%, 32 = AO 70%, 48 = AV 20%, 64 = AV 70%, 80 = VO

preStim = 0; % how much extra time do we want to include before trigger onset?
stimLength = 6; % onset still face between 1 and 1.3 seconds, offset still face is between 0.5 and 0.65 seconds,
% min video length = 5.54 s

cfg = [];

cfg.trialdef.eventvalue = triggerCode;
cfg.trialdef.prestim = preStim;
cfg.trialdef.poststim = stimLength;
cfg.channel    = {'MEGGRAD'};
cfg.continuous = 'yes';
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.eventtype = 'STI101'; % trigger channel name

for sub = 1:length(subjPaths),
    
    sprintf('=========================== Segmenting trials for %s ===========================\n', subjPaths{sub}(2:11))
    
    subjDir = strjoin(strcat(dataPath,subjPaths(sub)));
    cd(subjDir);
    
    cfg.dataset = strjoin(strcat(dataPath,subjPaths(sub),fileName_in)); % filename

    trials = ft_definetrial(cfg); % define trials

    trigger = cfg.trialdef.eventvalue;
        
    trials_allcond_grad = ft_preprocessing(trials); % get data based on the defined trials

    save trials_icaed_allcond_MEGGRAD.mat trials_allcond_grad
    
    % plot trial 5 on channel 100
    figure(sub);
    plot(trials_allcond_grad.time{5}, trials_allcond_grad.trial{5}(100,:)) % .time = seconds, .trial = data
    
end



% % plot trial 1 on channel 100
% figure(1);
% plot(trials_allcond.time{1}, trials_allcond.trial{1}(100,:)) % .time = seconds, .trial = data



