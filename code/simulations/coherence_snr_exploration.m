% Generating two sets of signals coherent with each-other
% Generate a signal consisting of a sinusoid signal. 
clearvars;

fs = 256;                    % Sampling frequency (samples per second)
dt = 1/fs;                   % seconds per sample
StopTime = 5;                % length of signal in seconds
t = (0:dt:StopTime-dt)';     % seconds
F = 4;                       % Frequency of sinusoid (Hz)
nSig = 100; 
nReps = 100;

% Generate 100 sinusoid signals with base frequency F and randomly set 
% amplitude and phase
ax = rand(nSig,nReps);
phx = 0.15*pi*randn(nSig,nReps);
sx = arrayfun(@(a,p) a*sin((2*pi*F*t)+p),ax,phx,'UniformOutput',false);
% Generate 100 sinusoid signals with base frequency F and randomly set 
% amplitude and phase
ay = rand(nSig,nReps);
phy = 0.15*pi*randn(nSig,nReps);
sy = arrayfun(@(a,p) a*sin((2*pi*F*t)+p),ay,phy,'UniformOutput',false);
% Add Gaussian noise to the sinusoids
snr_x = 1;
snr_y = snr_x;
x = cellfun(@(c) c + snr_x^-1*randn(fs*StopTime,1),sx,'UniformOutput',false);
x = arrayfun(@(i) cat(2,x{:,i}),1:nReps,'UniformOutput',false);
y = cellfun(@(c) c + snr_y^-1*randn(fs*StopTime,1),sy,'UniformOutput',false);
y = arrayfun(@(i) cat(2,y{:,i}),1:nReps,'UniformOutput',false);

% Define system responses
syst_sum = @(x,y) x + y;
syst_diff = @(x,y) x - y;

% Estimate the magnitude-squared coherence
[cxy,fc] = cellfun(@(a,b) coherence(a,b,fs),x,y,'UniformOutput',false);
fc = fc{1};
c_sumxy_x = cellfun(@(a,b) coherence(syst_sum(a,b),a,fs),x,y,'UniformOutput',false);
c_sumxy_y = cellfun(@(a,b) coherence(syst_sum(a,b),b,fs),x,y,'UniformOutput',false);
c_diffxy_x = cellfun(@(a,b) coherence(syst_diff(a,b),a,fs),x,y,'UniformOutput',false);
c_diffxy_y = cellfun(@(a,b) coherence(syst_diff(a,b),b,fs),x,y,'UniformOutput',false);
c_diffyx_x = cellfun(@(a,b) coherence(syst_diff(b,a),a,fs),x,y,'UniformOutput',false);
c_diffyx_y = cellfun(@(a,b) coherence(syst_diff(b,a),b,fs),x,y,'UniformOutput',false);

% Handy functions to compute mean and std across repetitions
fun_mean = @(x) mean(cat(2,x{:}),2);
fun_std = @(x) std(cat(2,x{:}),[],2);

cxy_mean = fun_mean(cxy);
cxy_std = fun_std(cxy);
c_sumxy_x_mean = fun_mean(c_sumxy_x);
c_sumxy_x_std = fun_std(c_sumxy_x);
c_sumxy_y_mean = fun_mean(c_sumxy_y);
c_sumxy_y_std = fun_std(c_sumxy_y);
c_diffxy_x_mean = fun_mean(c_diffxy_x);
c_diffxy_x_std = fun_std(c_diffxy_x);
c_diffxy_y_mean = fun_mean(c_diffxy_y);
c_diffxy_y_std = fun_std(c_diffxy_y);
c_diffyx_x_mean = fun_mean(c_diffyx_x);
c_diffyx_x_std = fun_std(c_diffyx_x);
c_diffyx_y_mean = fun_mean(c_diffyx_y);
c_diffyx_y_std = fun_std(c_diffyx_y);

% Plot results
figure();
set(gcf, 'Position', [50 50 1050 850], 'PaperPositionMode', 'auto');
subplot(3,4,1:4);
plot(t,x{1}(:,1)); hold on
plot(t,y{1}(:,1));
xlabel('time (s)');
ylabel('signal ampl [a.u.]');
legend('signal x1', 'signal y1');
title(sprintf('Example 4 Hz signals with Gaussian noise added, SNR:%.2f',snr_x));

subplot(3,4,6);
[u,v] = pol2cart(phx(:,1),ax(:,1));
compass(u,v);
title('Phase angles (x)');
subplot(3,4,7);
[u,v] = pol2cart(phy(:,1),ay(:,1));
compass(u,v);
title('Phase angles (y)');

subplot(3,4,9:10);
errorbar(fc,cxy_mean,cxy_std); hold on;
errorbar(fc,c_sumxy_x_mean,c_sumxy_x_std);
errorbar(fc,c_sumxy_y_mean,c_sumxy_y_std);
errorbar(fc,c_diffxy_x_mean,c_diffxy_x_std);
errorbar(fc,c_diffxy_y_mean,c_diffxy_y_std);
errorbar(fc,c_diffyx_x_mean,c_diffyx_x_std);
errorbar(fc,c_diffyx_y_mean,c_diffyx_y_std);
xlim([0,F+6]);
ylim([0,1]);
xlabel('frequency (Hz)');
ylabel('coherence');
title(sprintf('Coherence spectrum'));
legendStr = {'coh(x,y)','coh(x,x+y)','coh(y,x+y)','coh(x,x-y)','coh(y,x-y)',...
             'coh(x,y-x)','coh(y,y-x)'};
legend(legendStr);

subplot(3,4,11:12);
hzIdx = fc == F; 
barwitherr([cxy_std(hzIdx),c_sumxy_x_std(hzIdx),c_sumxy_y_std(hzIdx),...
            c_diffxy_x_std(hzIdx),c_diffxy_y_std(hzIdx),c_diffyx_x_std(hzIdx),...
            c_diffyx_y_std(hzIdx)],...
           [cxy_mean(hzIdx),c_sumxy_x_mean(hzIdx),c_sumxy_y_mean(hzIdx),...
            c_diffxy_x_mean(hzIdx),c_diffxy_y_mean(hzIdx),c_diffyx_x_mean(hzIdx),...
            c_diffyx_y_mean(hzIdx)]);
ylim([0,1]);
ylabel('coherence');
set(gca,'XTickLabel',legendStr);
xtickangle(45);
% legend({'sys(same)','sys(same)|p(other)','sys(same)|p(rand)','sys(other)',...
%         'sys(other)|p(other)'},'Location','southeast');
title(sprintf('Coherence at %d Hz',F));

%% 
% Generating two sets of signals coherent with each-other
% Generate a signal consisting of a sinusoid signal.
fs = 256;                    % Sampling frequency (samples per second)
dt = 1/fs;                   % seconds per sample
StopTime = 5;                % length of signal in seconds
t = (0:dt:StopTime-dt)';     % seconds
F = 4;                       % Frequency of sinusoid (Hz)
nSig = 100; 
nReps = 100;

% Generate 100 sinusoid signals with base frequency F and randomly set 
% amplitude and phase
ax = rand(nSig,1);
phx = 0.15*pi*randn(nSig,1);
sx = arrayfun(@(a,p) a*sin((2*pi*F*t)+p),ax,phx,'UniformOutput',false);
% Generate 100 sinusoid signals with base frequency F and randomly set 
% amplitude and phase
ay = rand(nSig,1);
phy = 0.15*pi*randn(nSig,1);
sy = arrayfun(@(a,p) a*sin((2*pi*F*t)+p),ay,phy,'UniformOutput',false);
% Add Gaussian noise to the sinusoids
snr_x = 1;
snr_y = snr_x;
x = cellfun(@(c) c + snr_x^-1*randn(fs*StopTime,1),sx,'UniformOutput',false);
x = cat(2,x{:});
y = cellfun(@(c) c + snr_y^-1*randn(fs*StopTime,1),sy,'UniformOutput',false);
y = cat(2,y{:});

% Define system responses
snr_s = 0.01;  % signal-to-noise ratio of the system
syst_id = @(x) x + snr_s^-1*randn(size(x));
syst_sum = @(x,y) x + y + snr_s^-1*randn(size(x));
syst_diff = @(x,y) x - y + snr_s^-1*randn(size(x));
% Compute system responses to input
sr_id_x = arrayfun(@(i) syst_id(x),1:nReps,'UniformOutput',false);
sr_id_y = arrayfun(@(i) syst_id(y),1:nReps,'UniformOutput',false);
sr_sum = arrayfun(@(i) syst_sum(x,y),1:nReps,'UniformOutput',false);
sr_diffxy = arrayfun(@(i) syst_diff(x,y),1:nReps,'UniformOutput',false);
sr_diffyx = arrayfun(@(i) syst_diff(y,x),1:nReps,'UniformOutput',false);
% Estimate the magnitude-squared coherence
[c_idx_x,fc] = cellfun(@(a) coherence(x,a,fs),sr_id_x,'UniformOutput',false);
fc = fc{1};
c_idy_y = cellfun(@(a) coherence(y,a,fs),sr_id_y,'UniformOutput',false);

c_sumxy_x = cellfun(@(a) coherence(x,a,fs),sr_sum,'UniformOutput',false);
c_sumxy_y = cellfun(@(a) coherence(y,a,fs),sr_sum,'UniformOutput',false);

c_diffxy_x = cellfun(@(a) coherence(x,a,fs),sr_diffxy,'UniformOutput',false);
c_diffxy_y = cellfun(@(a) coherence(y,a,fs),sr_diffxy,'UniformOutput',false);
c_diffyx_x = cellfun(@(a) coherence(x,a,fs),sr_diffyx,'UniformOutput',false);
c_diffyx_y = cellfun(@(a) coherence(y,a,fs),sr_diffyx,'UniformOutput',false);

% Handy functions to compute mean and std across repetitions
fun_mean = @(x) mean(cat(2,x{:}),2);
fun_std = @(x) std(cat(2,x{:}),[],2);

c_xx_mean = fun_mean(c_idx_x);
c_xx_std = fun_std(c_idx_x);
c_yy_mean = fun_mean(c_idy_y);
c_yy_std = fun_std(c_idy_y);
c_sumxy_x_mean = fun_mean(c_sumxy_x);
c_sumxy_x_std = fun_std(c_sumxy_x);
c_sumxy_y_mean = fun_mean(c_sumxy_y);
c_sumxy_y_std = fun_std(c_sumxy_y);
c_diffxy_x_mean = fun_mean(c_diffxy_x);
c_diffxy_x_std = fun_std(c_diffxy_x);
c_diffxy_y_mean = fun_mean(c_diffxy_y);
c_diffxy_y_std = fun_std(c_diffxy_y);
c_diffyx_x_mean = fun_mean(c_diffyx_x);
c_diffyx_x_std = fun_std(c_diffyx_x);
c_diffyx_y_mean = fun_mean(c_diffyx_y);
c_diffyx_y_std = fun_std(c_diffyx_y);

% Partial coherence
c_xx_py = cellfun(@(a) coherence(x,a,fs,'partialize',y),sr_id_x,'UniformOutput',false);
c_yy_px = cellfun(@(a) coherence(y,a,fs,'partialize',x),sr_id_y,'UniformOutput',false);
c_sumxy_x_py = cellfun(@(a) coherence(x,a,fs,'partialize',y),sr_sum,'UniformOutput',false);
c_sumxy_y_px = cellfun(@(a) coherence(y,a,fs,'partialize',x),sr_sum,'UniformOutput',false);
c_diffxy_x_py = cellfun(@(a) coherence(x,a,fs,'partialize',y),sr_diffxy,'UniformOutput',false);
c_diffxy_y_px = cellfun(@(a) coherence(y,a,fs,'partialize',x),sr_diffxy,'UniformOutput',false);

c_xx_py_mean = fun_mean(c_xx_py);
c_xx_py_std = fun_std(c_xx_py);
c_yy_px_mean = fun_mean(c_yy_px);
c_yy_px_std = fun_std(c_yy_px);
c_sumxy_x_py_mean = fun_mean(c_sumxy_x_py);
c_sumxy_x_py_std = fun_std(c_sumxy_x_py);
c_sumxy_y_px_mean = fun_mean(c_sumxy_y_px);
c_sumxy_y_px_std = fun_std(c_sumxy_y_px);
c_diffxy_x_py_mean = fun_mean(c_diffxy_x_py);
c_diffxy_x_py_std = fun_std(c_diffxy_x_py);
c_diffxy_y_px_mean = fun_mean(c_diffxy_y_px);
c_diffxy_y_px_std = fun_std(c_diffxy_y_px);

% Plot results
figure();
set(gcf, 'Position', [50 50 1050 850], 'PaperPositionMode', 'auto');
subplot(3,4,1:4);
plot(t,x(:,1)); hold on
plot(t,y(:,1));
xlabel('time (s)');
ylabel('signal ampl [a.u.]');
legend('signal x1', 'signal y1');
title(sprintf('Example 4 Hz signals with Gaussian noise added, SNR:%.2f',snr_x));

subplot(3,4,6);
[u,v] = pol2cart(phx,ax);
compass(u,v);
title('Phase angles (x)');
subplot(3,4,7);
[u,v] = pol2cart(phy,ay);
compass(u,v);
title('Phase angles (y)');

subplot(3,4,9:12);
hzIdx = fc == F; 
barwitherr([c_xx_std(hzIdx),c_xx_py_std(hzIdx),c_yy_std(hzIdx),c_yy_px_std(hzIdx),...
            c_sumxy_x_std(hzIdx),c_sumxy_x_py_std(hzIdx),...
            c_sumxy_y_std(hzIdx),c_sumxy_y_px_std(hzIdx)...
            c_diffxy_x_std(hzIdx),c_diffxy_x_py_std(hzIdx),...
            c_diffxy_y_std(hzIdx),c_diffxy_y_px_std(hzIdx)],...
           [c_xx_mean(hzIdx),c_xx_py_mean(hzIdx),c_yy_mean(hzIdx),c_yy_px_mean(hzIdx),...
            c_sumxy_x_mean(hzIdx),c_sumxy_x_py_mean(hzIdx),...
            c_sumxy_y_mean(hzIdx),c_sumxy_y_px_mean(hzIdx),...
            c_diffxy_x_mean(hzIdx),c_diffxy_x_py_mean(hzIdx),...
            c_diffxy_y_mean(hzIdx),c_diffxy_y_px_mean(hzIdx)]);
% ylim([0,1]);
ylabel('coherence');
legendStr = {'coh(x,syst(x))','coh(x,syst(x)|py)','coh(y,syst(y))',...
             'coh(y,syst(y)|px)','coh(x,syst(x+y))','coh(x,syst(x+y)|py)',...
             'coh(y,syst(x+y))','coh(y,syst(x+y)|px)','coh(x,syst(x-y))',...
             'coh(x,syst(x-y)|py)','coh(y,syst(x-y))','coh(y,syst(x-y)|px)'};
set(gca,'XTick',1:numel(legendStr),'XTickLabel',legendStr);
xtickangle(45);
% legend({'sys(same)','sys(same)|p(other)','sys(same)|p(rand)','sys(other)',...
%         'sys(other)|p(other)'},'Location','southeast');
title(sprintf('Coherence at %d Hz, system SNR: %.2f',F,snr_s));
