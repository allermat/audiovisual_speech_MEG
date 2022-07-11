%% Compute coherence between auditory and visual signals
df_path = fullfile(AVSM_setupdir('analysis_megcoherence_group'),...
                                 'ftmeg_coh_acoustic_lip_truncated.mat');
% Only compute if the results are not ready
if ~exist(df_path,'file')
    coh_acoustic_lip();
end

%% Plot figure
% Load the data
load(fullfile(AVSM_setupdir('analysis_megcoherence_group'),...
     'ftmeg_coh_acoustic_lip_truncated.mat'),'coh','cohPerm','ftData_pow',...
     'ftData_pow_baseline');

% sort the permuted coherences for each frequency bin separately
cohPermSorted = sort(cohPerm,2);

% find the boudaries of 95%, 99.5%, 99.9% CI
idx95 = (size(cohPermSorted,2)*0.95);
idx995 = (size(cohPermSorted,2)*0.995);
idx999 = (size(cohPermSorted,2)*0.999);
idx99875 = round(size(cohPermSorted,2)*0.99875); % Bonferroni corrected (40 tests) p < 0.05
ci95 = cohPermSorted(:,idx95);
ci995 = cohPermSorted(:,idx995);
ci999 = cohPermSorted(:,idx999);
ci99875 = cohPermSorted(:,idx99875); % Bonferroni corrected (40 tests) p < 0.05

% Plot powerspectra
figure('Position',[100 100 1400 400]);
subplot(1,3,1);
plot(ftData_pow.freq,squeeze(mean(ftData_pow.powspctrm)),'LineWidth',1.5); hold on;
plot(ftData_pow_baseline.freq,squeeze(mean(ftData_pow_baseline.powspctrm)),'LineWidth',1.5); hold on;
legend({'Auditory','Lip','Auditory baseline','Lip baseline'});
title(sprintf('Auditory and Lip envelope powerspectra\ntruncated sentences (5s)'));
xlabel('Frequency (Hz)');
ylabel('Power');

% log-log plot of powerspectra
subplot(1,3,2);
% The fitting is based on Chandrasekaran et al. 2009
% The fitted function is Y=a*f^b in the logarithmic form
% log(Y)=log(a)-b*log(f). This way linear methods can be used for fitting. 
fitOneOverF = @(x,y) fit(x,y,fittype('a-b*x'));
logFreq = log10(ftData_pow.freq');
logSpctrm = log10(transpose(squeeze(mean(ftData_pow.powspctrm))));
audFit = fitOneOverF(logFreq,logSpctrm(:,1));
lipFit = fitOneOverF(logFreq,logSpctrm(:,2));
xx = linspace(min(logFreq),max(logFreq),200);

plot(logFreq,logSpctrm,'LineWidth',1.5); hold on
temp = cellstr(num2str(ftData_pow.freq(2:2:end)'));
xtickstr = repmat({''},size(temp));
xtickstr([1,2,4,8,16,20]) = temp([1,2,4,8,16,20]);
set(gca,'XTick',logFreq(2:2:end),'XTickLabel',xtickstr);
plot(xx,audFit(xx),'LineWidth',1.5,'LineStyle','--','Color',[0,0.4470,0.7410]);
plot(xx,lipFit(xx),'LineWidth',1.5,'LineStyle','--','Color',[0.8500,0.3250,0.0980]);

legend({'Auditory','Lip','1/f fit Aud','1/f fit Lip'});
title(sprintf('Auditory and Lip envelope powerspectra\ntruncated sentences (5s)'));
xlabel('Frequency (Hz)');
ylabel('Power (log10 units)');

% Plot coherence spectrum
subplot(1,3,3);
plot(coh.freq,squeeze(coh.cohspctrm(2,1,:)),'LineWidth',1.5); hold on;
plot(coh.freq,ci99875,'--r','LineWidth',1);
% plot(coh.freq,ci995,'--g','LineWidth',1);
% plot(coh.freq,ci999,'--c','LineWidth',1);
legend({'aud-lip coherence','p < 0.05 Bonferroni-corrected'});
title(sprintf('Auditory-Lip Coherence\ntruncated sentences (5s)'));
xlabel('Frequency (Hz)');
ylabel('Coherence');

