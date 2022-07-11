%% Visual artifact rejection
%
% Heidi Solberg Okland
%
% Use FieldTrip function ft_rejectvisual to reject trials that have either
% high variance or have amplitudes that are 4SD above/below the mean.

clear all; close all;

addpath(genpath('/imaging/local/software/fieldtrip/fieldtrip-latest'));

dataPath = '';

subjPaths = {};

triggerCode = [16 32 48 64 80];
% key: 16 = AO 20%, 32 = AO 70%, 48 = AV 20%, 64 = AV 70%, 80 = VO

Cond = {'AO 20%', 'AO 70%', 'AV 20%', 'AV 70%', 'VO'};

rawFile = 'trials_icaed_allcond_MEG';

% define outlier cutoff in standard deviation
cutoff = 4; 

percent_rejected_trials = {};

for sub = 1:length(subjPaths),
    
    current_sub = subjPaths{sub}(2:11);
    fprintf('=========================== Now looking at %s ===========================\n', current_sub)
    
    subjDir = strjoin(strcat(dataPath,subjPaths(sub)));
    cd(subjDir);
    
    % load data
    load(rawFile)
    
    % --------Reject trials based on extreme amplitudes--------
    
    % first calculate averaged trials (ERFs)
    cfg               = [];
    cfg.vartrllength  = 0;
    cfg.keeptrials    = 'yes';
    erf_trials = ft_timelockanalysis(cfg,trials_allcond);
    
    % calculate the average amplitude across time for each trial
    cfg = [];
    cfg.avgoverchan = 'yes';
    cfg.avgovertime = 'yes';
    avg_ampl = ft_selectdata(cfg, erf_trials);
    
    % calculate the amplitude that would be a certain number (cutoff) of 
    % standard deviations above the mean accross trials
    outlier_cutoff = cutoff*std(avg_ampl.trial);
    
    % find trials that have average amplitudes deviating more than the set
    % cutoff
    outlier_trials1 = find(avg_ampl.trial > outlier_cutoff);
    outlier_trials2 = find(avg_ampl.trial < -outlier_cutoff);
    outlier_trials = [outlier_trials1; outlier_trials2];
    
    outlier_trials % print outlier trials to the command line
    % Make sure you select the outlier trials to be rejected in the GUI!
    
    % --------Visually inspect trials to reject them--------
    
    % Define parameters
    cfg               = [];
    cfg.metric        = 'var';
    cfg.layout        = 'neuromag306all.lay';
    cfg.keepchannel   = 'yes'; % keeps channels that are not selected (e.g. keeps gradiometers when selecting to do magnetometers)
    cfg.keeptrial     = 'nan'; % replace bad trials with NaNs (so that we can see which ones were rejected)

    % Magnetometers
    cfg.channel       = 'MEGMAG';
    trials_clean_temp    = ft_rejectvisual(cfg,trials_allcond);
    
    % Gradiometers
    cfg.channel       = 'MEGGRAD';
    trials_clean_allcond    = ft_rejectvisual(cfg,trials_clean_temp);
    
    % --------Find out which trials were rejected---------------
    
    nTrials = length(trials_clean_allcond.trial);
    cond = trials_clean_allcond.trialinfo;
    trials_kept = [];
    trials_rejected = [];
    keepTrial = nan(2,nTrials);
    
    for tr = 1:nTrials
    
        % Pull out the relevant trial
        thisTrial = trials_clean_allcond.trial{1,tr};
            
        % Check if it's got NaNs in order to know which trials were
        % rejected
        if isnan(thisTrial)       
            fprintf('Trial %d is marked as BAD.\n', tr)
            trials_rejected = [trials_rejected; tr];
            keepTrial(1,tr) = 0; % bad trials marked with 0
        elseif ~isnan(thisTrial)
%             fprintf('Trial %d is marked as good.\n', tr)
            trials_kept = [trials_kept tr];
            keepTrial(1,tr) = 1; % good trials marked with 1
        end
        
        keepTrial(2,tr) = cond(tr); % Make a note of which condition the good and bad trials belong to
        
    end
    
    percent_rejected = (length(trials_rejected)/nTrials)*100;
    fprintf('\n### %.2f percent of trials rejected. ###\n\n', percent_rejected)
    
    percent_rejected_trials{sub,1} = current_sub;
    percent_rejected_trials{sub,2} = percent_rejected;
    
    keepTrial = keepTrial';
    
%     nRejected = length(trials_rejected);
%     oldTrials = trials_clean_allcond;
%     
%     % Remove trials that were marked as bad, otherwise just save
%     if nRejected > 0
%         
%     % Remove the rejected trials and info associated with them
%     fprintf('Deleting %d trials...\n', nRejected)
%     
%     CleanTrials = {};
%     SampleInfo = [];
%     Time = {};
%     TrialInfo = [];
%     for i = 1:length(trials_kept),
%         CleanTrials{1,i} = trials_clean_allcond.trial{1,trials_kept(i)};
%         SampleInfo(i,1) = trials_clean_allcond.sampleinfo(trials_kept(i),1);
%         SampleInfo(i,2) = trials_clean_allcond.sampleinfo(trials_kept(i),2);
%         Time{1,i} = trials_clean_allcond.time{1,trials_kept(i)};
%         TrialInfo(i) = trials_clean_allcond.trialinfo(trials_kept(i));
%     end
%     
%     % Remove rejected trials
%     trials_clean_allcond.trial = CleanTrials;
%     trials_clean_allcond.sampleinfo = SampleInfo;
%     trials_clean_allcond.time = Time;
%     trials_clean_allcond.trialinfo = TrialInfo';
%     
%     elseif nRejected == 0
%         fprintf('No bad trials found, so not rejecting any trials.')
%     end

    
    % Save!
    fprintf('Saving...\n')
    save trials_clean_allcond_MEG trials_clean_allcond
    save trials_rejected trials_rejected
    save trials_kept trials_kept
    save keepTrial keepTrial
    fprintf('\n=========================== Done with %s ===========================\n\n', subjPaths{sub}(2:11))
    
end

cd(dataPath)
save percent_rejected_trials percent_rejected_trials

% %% 3) EEG
% 
% cfg               = [];
% cfg.metric        = 'var';
% cfg.channel       = 'EEG';
% % cfg.layout        = ''; NEEDS SPECIFICATION!
% cfg.keepchannel   = 'yes';  % This keeps those channels that are not displayed in the data
% 
% for sub = 1:length(subjPaths),
%     
%     sprintf('======Now looking at EEG for subject %s======\n', subjPaths{1}(1:10))
%     subjDir = strjoin(strcat(dataPath,subjPaths(sub)));
%     cd(subjDir);
%     
%     for trig = 1:length(Cond), % for every condition
%         
%         trigger = triggerCode(trig);
%         fprintf('======Processing trials from the %s condition======\n', strjoin(Cond(trig)));
% 
%         if trigger == 16,
%             trials_AO20_clean2    = ft_rejectvisual(cfg,trials_AO20_clean);
%         elseif trigger == 32,
%             trials_AO70_clean2    = ft_rejectvisual(cfg,trials_AO70_clean);
%         elseif trigger == 48,
%             trials_AV20_clean2    = ft_rejectvisual(cfg,trials_AV20_clean);
%         elseif trigger == 64,
%             trials_AV70_clean2    = ft_rejectvisual(cfg,trials_AV70_clean);
%         elseif trigger == 80,
%             trials_VO_clean2    = ft_rejectvisual(cfg,trials_VO_clean);
%         end
%         
%     end
%     
% end

% % Save!
% save('trials_clean_MEEG.mat', 'trials_AO20_clean2', 'trials_AO70_clean2','trials_AV20_clean2', 'trials_AV70_clean2','trials_VO_clean2')
