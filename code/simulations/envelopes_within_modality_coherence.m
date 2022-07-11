clearvars;
%% Using the code from my simulations
% Load data
sourceDir = AVSM_setupdir('analysis_megcoherence_sub','sub-01');
data = load(fullfile(sourceDir,'trials_with_env_MEG'));
data = data.trials_with_env_MEG;
audIdx = ismember(data.label,'AudEnv');
visIdx = ismember(data.label,'LipEnv');
aud_env = cellfun(@(c) c(audIdx,:)',data.trial,'UniformOutput',false);
vis_env = cellfun(@(c) c(visIdx,:)',data.trial,'UniformOutput',false);

% Compute all possible pairings of sentences, leaving pairs between same
% sentences out
pairings = nchoosek(1:numel(aud_env),2);
x_aud = cat(2,aud_env{pairings(:,1)});
y_aud = cat(2,aud_env{pairings(:,2)});
x_vis = cat(2,vis_env{pairings(:,1)});
y_vis = cat(2,vis_env{pairings(:,2)});

[cxy_aud,fc] = coherence(x_aud,y_aud,data.fsample,'square',false);
cxy_vis = coherence(x_vis,y_vis,data.fsample,'square',false);

% Plot results
figure(); 
plot(fc,cxy_aud); hold on;
plot(fc,cxy_vis);
xlim([0,20]);
ylim([0,1]);
legend('Auditory','Visual');
xlabel('frequency (Hz)');
ylabel('coherence');
title('Within-modality coherence across sentences');

% Compue coherence spectrum between aud-vis for good measure
cxy_aud_vis = coherence(cat(2,aud_env{:}),cat(2,vis_env{:}),data.fsample,'square',false);

% Plot results
figure(); 
plot(fc,cxy_aud_vis);
xlim([0,20]);
% ylim([0,1]);
xlabel('frequency (Hz)');
ylabel('coherence');
title('Auditory-lip coherence');

%% Using fieldtrip
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

% Compute all possible pairings of sentences, leaving pairs between same
% sentences out
pairings = nchoosek(1:numel(aud_env),2);
% Feed all possible combinations sentences into data structure, auditory
% envelopes
ftData_f_allpairs_aud = ftData_fourier;
temp1 = cat(1,ftData_fourier.fourierspctrm(pairings(:,1),1,:));
temp2 = cat(1,ftData_fourier.fourierspctrm(pairings(:,2),1,:));
ftData_f_allpairs_aud.fourierspctrm = cat(2,temp1,temp2);
ftData_f_allpairs_aud.cumsumcnt = repmat(ftData_fourier.cumsumcnt(1,:),...
    size(pairings,1),1);
ftData_f_allpairs_aud.cumtapcnt = repmat(ftData_fourier.cumtapcnt(1,:),...
    size(pairings,1),1);
ftData_f_allpairs_aud = rmfield(ftData_f_allpairs_aud,'trialinfo');
cfg = [];
cfg.method = 'coh';
cfg.channelcmb = {'AudEnv' 'LipEnv'};
coh_aud = ft_connectivityanalysis(cfg, ftData_f_allpairs_aud);

% Feed all possible combinations sentences into data structure, auditory
% envelopes
ftData_f_allpairs_vis = ftData_fourier;
temp1 = cat(1,ftData_fourier.fourierspctrm(pairings(:,1),2,:));
temp2 = cat(1,ftData_fourier.fourierspctrm(pairings(:,2),2,:));
ftData_f_allpairs_vis.fourierspctrm = cat(2,temp1,temp2);
ftData_f_allpairs_vis.cumsumcnt = repmat(ftData_fourier.cumsumcnt(1,:),...
    size(pairings,1),1);
ftData_f_allpairs_vis.cumtapcnt = repmat(ftData_fourier.cumtapcnt(1,:),...
    size(pairings,1),1);
ftData_f_allpairs_vis = rmfield(ftData_f_allpairs_vis,'trialinfo');
cfg = [];
cfg.method = 'coh';
cfg.channelcmb = {'AudEnv' 'LipEnv'};
coh_vis = ft_connectivityanalysis(cfg, ftData_f_allpairs_vis);

figure();
plot(coh_aud.freq,squeeze(coh_aud.cohspctrm(2,1,:)),'LineWidth',1.5); hold on;
plot(coh_vis.freq,squeeze(coh_vis.cohspctrm(2,1,:)),'LineWidth',1.5);
ylim([0,1]);
legend({'Auditory','Visual'});
title(sprintf('Withon-modality coherence across sentences'));
xlabel('frequency (Hz)');
ylabel('coherence');