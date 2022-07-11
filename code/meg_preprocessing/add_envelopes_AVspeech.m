%% Trim trials and add corresponding speech envelopes
%
% Heidi Solberg Okland
%
% This script inserts speech envelopes and lip aperture timecourses into
% the relevant MEG trial based on the logfile made during the experiment.
% It also trims trials so that they are all the same length.
% Run this script after finishing preprocessing and cleaning the data.

clear all; close all; clc;

addpath(genpath('/imaging/local/software/fieldtrip/fieldtrip-latest'));

dataPath = '';

subjPaths = {};

inputData = 'trials_clean_allcond_GRAD';
rejectedTrials = 'trial_rejection_GRAD';
% Triggers: 16 = AO 20%, 32 = AO 70%, 48 = AV 20%, 64 = AV 70%, 80 = VO

load('BaselineGroupIndex.mat')
load('AudEnv_250Hz.mat')
load('LipEnv_250Hz.mat')

trial_length = 5;

cfg = [];

for sub = 1:length(subjPaths),
    
    % Start the timer
    tic 
    fprintf('=========================== Now adding envelopes for %s ===========================\n', subjPaths{sub}(2:11))
    
    % Go to subject directory
    subjDir = strjoin(strcat(dataPath,subjPaths(sub)));
    cd(subjDir);
    % Load logfile that contains info about which sentences were presented in which trial
    load logfile % (variable name: Trials)
    
    % Load MEG data and info about rejected trials
    load(inputData);
    load(rejectedTrials);
    oldData = trials_clean_allcond_grad;
    Fs = oldData.fsample;
    nTrialsOld = length(trials_clean_allcond_grad.trial);
    nTrialsNew = length(trials_kept);
    
    % -----------Make replacement variables-----------------
    
    newData = [];
    
    % Insert things we don't need to change
    newData.hdr = oldData.hdr;
    newData.fsample = oldData.fsample;
    newData.grad = oldData.grad;
    newData.elec = oldData.elec;
    newData.cfg = oldData.cfg;
    
    % Make new empty variables for things that need to change
    newData.label = {};
    newData.trial = {};
    newData.trialinfo = [];
    newData.time = {};
    newData.sampleinfo = [];
    
    % Make new, extra channel names for the envelopes, and insert into "label" field
   
    nMEGch = length(oldData.label);
    newLabels = cell(nMEGch+2,1);
    chAudEnv = nMEGch+1;
    chLipEnv = nMEGch+2;

    newLabels(1:nMEGch,1) = oldData.label;
    newLabels(chAudEnv,1) = {'AudEnv'};
    newLabels(chLipEnv,1) = {'LipEnv'};
    
    newData.label = newLabels; % insert!

        
        % Loop over trials to trim and insert envelopes + change associated
        % time vectors etc.

        newTrialCounter = 0;
        
            NaN_LipEnv = [];
            NaN_AudEnv = [];
            NaN_Trial = [];
            NaN_found = 0;
            
        for tr = 1:nTrialsOld,
            
            % Check if current trial is marked as bad ..
            if keepTrial(tr) == false,
                
               fprintf('*Trial %d is marked as bad. Skipping this one.*\n', tr)
               
            % Check if current trial is marked as good ..
            elseif keepTrial(tr) == true, 
                
            newTrialCounter = newTrialCounter+1;

            % Pull out the relevant trial
            oldTrial = oldData.trial{1,tr};
            
            % Check if it's got NaNs (just in case..)
            if isnan(oldTrial)       
                error('Trial %d is apparently marked as "good" but consists of NaNs..\nSomething is wrong :(\n', tr)
            end

            % Get the sentence index from the logfile
            sentenceID = char(Trials{tr,3}); % convert cell array to characters
            sentenceIdx = str2num(sentenceID(2:4)); %.. and convert characters to numbers!

            % Get the right envelopes based on this sentence index
            oldLipEnv = LipEnv_250Hz{1,sentenceIdx};
            oldAudioEnv = AudEnv_250Hz{1,sentenceIdx};

            % Get the pre- and poststim info for this sentence and calculate by
            % how much the trial needs to be shifted forwards
            
            if BaselineGroupIndex(sentenceIdx,2) == 1
            PreBaseline = 1.0;
            PostBaseline = 0.65;
            ShiftForward = 0;

            elseif BaselineGroupIndex(sentenceIdx,2) == 2 
            PreBaseline = 1.1;
            PostBaseline = 0.6;
            ShiftForward = PreBaseline-1;

            elseif BaselineGroupIndex(sentenceIdx,2) == 3
            PreBaseline = 1.2;
            PostBaseline = 0.55;
            ShiftForward = PreBaseline-1;

            elseif BaselineGroupIndex(sentenceIdx,2) == 4
            PreBaseline = 1.3;
            PostBaseline = 0.5;
            ShiftForward = PreBaseline-1;

            end

                % Shift and trim trials & envelopes

                if ShiftForward > 0, % if we need to shift things..

                % Find out where to cut..
                firstSample = round(ShiftForward*Fs)+1;
                lastSample = trial_length*Fs+firstSample-1;

                % Shift the trial onset and trim trials and envelopes to desired length
                newTrial = oldTrial(:,firstSample:lastSample);
                newLipEnv = oldLipEnv(firstSample:lastSample);
                newAudEnv = oldAudioEnv(firstSample:lastSample);

                elseif ShiftForward == 0, % if we don't..

                % Trim the the trial and envelopes to the desired length
                lastSample = trial_length*Fs;
                newTrial = oldTrial(:,1:lastSample);
                newLipEnv = oldLipEnv(1:lastSample);
                newAudEnv = oldAudioEnv(1:lastSample);

                end
        
            % Check that you've done the math right..
            oldLength = length(oldTrial);
            newLength = length(newTrial);
            if newLength ~= trial_length*Fs
               error('Trial %d has now been trimmed such that it is now %d samples long, but should have been %d samples.', tr, newLength, trial_length*Fs)
            end
            if oldLength == newLength
                warning('Trial %d is still exactly the same length. Looks like something has gone wrong...', tr)
            end
          
            % Check if the envelopes or trials have NaNs - if so, give error!
            if find(isnan(newTrial)) > 0
                warning('Trial %d seems to have NaNs in it after trimming. Something is wrong :(', newTrialCounter)
                NaN_Trial = [NaN_Trial ; newTrialCounter];
                NaN_found = NaN_found+1;
            end
          
            if find(isnan(newLipEnv)) > 0
                
                warning('The lip aperture for new trial %d has NaNs in it after trimming. Something is wrong :(',newTrialCounter)
                NaN_LipEnv = [NaN_LipEnv ; newTrialCounter];
                NaN_found = NaN_found+1;
                
                figure;
                MEGsens = 50;
                fprintf('Plotting MEG sensor %d and envelopes for sentence %s...\n', MEGsens, sentenceID)
                plot(newData.time{1,1}(1,:), newTrial(MEGsens,:), 'y'); hold on;
                plot(newData.time{1,1}(1,:), newAudEnv*10^-10, 'r'); hold on;
                plot(newData.time{1,1}(1,:), newLipEnv*10^-15, 'b')
                title('Trial with NaNs in lip aperture timecourse')
                
            end
            
            if find(isnan(newAudEnv)) > 0
                warning('The audio envelope for new trial %d has NaNs in it after trimming. Something is wrong :(',newTrialCounter)
                NaN_AudEnv = [NaN_AudEnv ; newTrialCounter];
                NaN_found = NaN_found+1;
            end
            
            
            % Insert new trial
            expandedTrial = zeros(nMEGch+2,newLength);
            expandedTrial(1:nMEGch,:) = newTrial;
            expandedTrial(chAudEnv,:) = newAudEnv*10^-10; % scale the audio envelope down
            expandedTrial(chLipEnv,:) = newLipEnv*10^-15; % scale the lip envelope down
            newData.trial{1,newTrialCounter} = expandedTrial; % insert!
            
            % Add condition index
            newData.trialinfo(newTrialCounter,1) = oldData.trialinfo(tr,1);
            
            % Change the "time" field
            newTime = oldData.time{1,1}(1:trial_length*Fs);
            newData.time{newTrialCounter} = newTime; % insert!
            
            % Change the "sampleinfo" field
            samplesShifted = round(ShiftForward * Fs);
            samplesTrimmed = oldLength-newLength;   

            oldStart = oldData.sampleinfo(tr,1);
            oldFinish = oldData.sampleinfo(tr,2);
            newStart = oldStart + samplesShifted;
            newFinish = oldFinish - samplesTrimmed + samplesShifted;

                % Again check the math...
                if newFinish - newStart + 1 ~= newLength
                error('The new sampleinfo does not match up with the length of the trimmed trial.')
                end
                
            % .. and insert!
            newData.sampleinfo(newTrialCounter,1) = newStart;
            newData.sampleinfo(newTrialCounter,2) = newFinish;
            
            end
            
        end
        
    % Give error if we found NaNs
    if NaN_found > 0,
        
    error('NaNs found in some trials, lip or audio envelopes. See variables starting with "NaN_" for trial indices.')
    
    end

    % Check that we don't have too few or too many trials now..
    if length(newData.trial) > nTrialsNew || length(newData.trial) < nTrialsNew,
        error('You now have %d trials in your struct, but you should have had %d.\nSomething is wrong!\n', length(newData.trial), nTrialsNew)
    end
    
    % Change the variable name
    trials_with_env_grad = newData;
   
    % Save in subject directory
    fprintf('Saving...\n')    
    save trials_with_env_MEGGRAD trials_with_env_grad
    
    % Stop the timer
    toc
    fprintf('=========================== Done with %s ===========================\n', subjPaths{sub}(2:11))
    
end



%% Plot things

% close all;
% 
% i = 1;
% for trial = 150:160,
% figure(i), hold on;
% plot(1:newLength,trials_with_env_MEG.trial{1,trial}(306,:), 'k');
% plot(1:newLength,trials_with_env_MEG.trial{1,trial}(307,:), 'r');
% plot(1:newLength,trials_with_env_MEG.trial{1,trial}(308,:), 'g');
% i = i+1;
% end

