%% Resample envelopes for audio and lip movements
%
% Heidi Solberg Økland, May 2017
%
% This script uses the function 'resample' to up/downsample speech
% envelopes and lip aperture timecourses. It pads the timecourses with 1s
% on each side of the signal in order to avoid artifacts.


clear all; close all;

% Load data matrices
load LipEnvelopes_MEG.mat
load AudioEnvelopes_MEG.mat

% Define sampling frequencies
Fs_vid = 50;
Fs_aud = 48000;
Fs_target = 250;
nSentences = 301;

amp = 10^5; % amplification factor for audio envelopes so as to compare with lip aperture

% Make empty cell arrays for the two data types
LipEnv_250Hz = cell(1,nSentences);
AudEnv_250Hz = cell(1,nSentences);

% Loop over sentences
for sent = 1:nSentences,
    
    % Get data for this sentence and remove NaNs
    AUD = AudioEnvelopes_MEG(sent,:);
    AUD(isnan(AUD))=[];
    LIP = LipEnvelopes_MEG(sent,:);
    LIP(isnan(LIP))=[];
    
    % Pad data with one second on each side
    AUDpad = [zeros(1,Fs_aud) AUD zeros(1,Fs_aud)];
    LIPpad = [repmat(LIP(1),1,Fs_vid) LIP repmat(LIP(end),1,Fs_vid)];
    
    % Resample!
    AUD250 = resample(AUDpad',Fs_target,Fs_aud);
    LIP250 = resample(LIPpad',Fs_target,Fs_vid);
    
    % Transpose
    AUD250 = AUD250';
    LIP250 = LIP250';
    
    % Remove padding
    AUD250 = AUD250(Fs_target+1:end-Fs_target);
    LIP250 = LIP250(Fs_target+1:end-Fs_target);
    
    % Insert into cell array
    LipEnv_250Hz{1,sent} = LIP250;
    AudEnv_250Hz{1,sent} = AUD250;
    
    % Check that they're the same length
    if length(LipEnv_250Hz{1,sent}) ~= length(AudEnv_250Hz{1,sent})
        warning('The lengths of the audio envelope and lip aperture timecourse for sentence %d do not match. Something is wrong... :(', sent)
        fprintf('Plotting sentence %d...\n', sent)
        figure; plot(AudEnv_250Hz{1,sent}*amp, 'r'), hold on; plot(LipEnv_250Hz{1,sent}, 'b')
    end
    
        % Plot the first sentence and check whether we want to move on
        if sent == 1,

        figure; plot(AudEnv_250Hz{1,sent}*amp, 'r'), hold on; plot(LipEnv_250Hz{1,sent}, 'b')

        prompt = 'OK to continue? (y/n)\n';
        answer = input(prompt,'s');

        if answer == 'y',
            fprintf('OK, great, continuing with the rest of the sentences...\n')
        elseif answer == 'n',
            fprintf('Not continuing.\n')
            return
        end
        end
    
end


% Happy? Then let's save!

prompt = 'OK to save? (y/n)\n';
answer = input(prompt,'s');

if answer == 'y',
    
fprintf('Saving...\n')
cd /imaging/hs01/AV_speech/
save AudEnv_250Hz AudEnv_250Hz
save LipEnv_250Hz LipEnv_250Hz

elseif answer == 'n',
    
fprintf('Not saving.\n')

end

