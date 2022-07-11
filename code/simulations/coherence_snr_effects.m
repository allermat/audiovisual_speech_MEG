% Simulations to determine stimulus and signal SNR to key coherence
% results
clearvars;
currPool = gcp('nocreate');
if isempty(currPool)
    parpool('local');
end
% Preparing the parameter space
% grid resolution
noise_stim = [0,0.25,0.5,1,2,4,8];
noise_syst = [0,0.25,0.5,1,2,4,8];
% Parameter names
parameterNames = {'noise_stim','noise_syst'};
gridVectors = {noise_stim,noise_syst};
nParameters = numel(gridVectors);
% Full factorial expansion of the specified parameter vectors
coords = cell(1,nParameters);
[coords{:}] = ndgrid(gridVectors{:});
coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
% Matrix of all possible parameter combinations: nRows = number of
% combinations, columns correspond to parameters
paramCombinations = cat(2,coords{:});

% Generating two sets of signals coherent with each-other
% Generate a signal consisting of a sinusoid signal.
fs = 256;                    % Sampling frequency (samples per second)
dt = 1/fs;                   % seconds per sample
StopTime = 5;                % length of signal in seconds
t = (0:dt:StopTime-dt)';     % seconds
F = 4;                       % Frequency of sinusoid (Hz)
nSig = 100; 
nSubj = 100;

% Handy functions to compute mean and std across repetitions
fun_mean = @(x) mean(cat(2,x{:}),2);
fun_std = @(x) std(cat(2,x{:}),[],2);

% Array for storing outputs
coh_all = cell(size(paramCombinations,1),4);
parfor_progress(size(paramCombinations,1));
parfor iGrid = 1:size(paramCombinations,1)
    
    noise_x = paramCombinations(iGrid,1); % Amount of noise added to the signals
    noise_s = paramCombinations(iGrid,2); % Amount of noise added by the system
       
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
    noise_y = noise_x;
    x = cellfun(@(c) c + noise_x*randn(fs*StopTime,1),sx,'UniformOutput',false);
    x = cat(2,x{:});
    y = cellfun(@(c) c + noise_y*randn(fs*StopTime,1),sy,'UniformOutput',false);
    y = cat(2,y{:});
    
    % Define system responses
    syst_id = @(x) x + noise_s*randn(size(x));
    syst_sum = @(x,y) x + y + noise_s*randn(size(x));
    
    % Compute system responses to input
    sr_id_x = arrayfun(@(i) syst_id(x),1:nSubj,'UniformOutput',false);
    sr_sum = arrayfun(@(i) syst_sum(x,y),1:nSubj,'UniformOutput',false);
    % Estimate the magnitude-squared coherence
    [c_idx_x,fc] = cellfun(@(a) coherence(x,a,fs),sr_id_x,'UniformOutput',false);
    fc = fc{1};
    c_xx_py = cellfun(@(a) coherence(x,a,fs,'partialize',y),sr_id_x,'UniformOutput',false);
    c_sumxy_x = cellfun(@(a) coherence(x,a,fs),sr_sum,'UniformOutput',false);
    c_sumxy_x_py = cellfun(@(a) coherence(x,a,fs,'partialize',y),sr_sum,'UniformOutput',false);

    hzIdx = fc == F;
    temp = {cat(2,c_idx_x{:})',cat(2,c_xx_py{:})',cat(2,c_sumxy_x{:})',...
            cat(2,c_sumxy_x_py{:})'};
    coh_all(iGrid,:) = cellfun(@(x) x(:,hzIdx),temp,'UniformOutput',false);
    
    parfor_progress;
end
parfor_progress(0);

% Plot figures
coh_all_tbl = cell2table(coh_all,'VariableNames',{'x_sIdx','x_sIdx_py',...
                         'x_sSumxy','x_sSumxy_py'});
X = reshape(paramCombinations(:,1),numel(noise_stim),numel(noise_syst));
Y = reshape(paramCombinations(:,2),numel(noise_stim),numel(noise_syst));
temp = cellfun(@minus,coh_all_tbl.x_sIdx,coh_all_tbl.x_sIdx_py,...
               'UniformOutput',false);
Z = reshape(cellfun(@mean,temp),numel(noise_stim),numel(noise_syst));
E = reshape(cellfun(@std,temp),numel(noise_stim),numel(noise_syst));
figure();
subplot(1,2,1);
surf(X,Y,Z); hold on
plot3([X(:),X(:)]', [Y(:),Y(:)]', [-E(:),E(:)]'+Z(:)','Color','k',...
      'LineWidth',1.5);
hold off
legend('mean','SD')
xlabel('Stimulus noise (SD)');
ylabel('System noise (SD)');
zlabel('Coherence difference');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
title('coh(x,sysID(x)) > coh(x,sysID(x)|part y)');
subplot(1,2,2);
contour(X,Y,Z);
xlabel('Stimulus noise (SD)');
ylabel('System noise (SD)');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
c = colorbar;
c.Label.String = 'Coherence difference';
colormap(crameri('vik'));


temp = cellfun(@minus,coh_all_tbl.x_sSumxy,coh_all_tbl.x_sIdx,...
               'UniformOutput',false);
Z = reshape(cellfun(@mean,temp),numel(noise_stim),numel(noise_syst));
E = reshape(cellfun(@std,temp),numel(noise_stim),numel(noise_syst));
figure();
subplot(1,2,1);
surf(X,Y,Z); hold on
plot3([X(:),X(:)]', [Y(:),Y(:)]', [-E(:),E(:)]'+Z(:)','Color','k',...
      'LineWidth',1.5);
hold off
legend('mean','SD')
xlabel('Stimulus noise (SD)');
ylabel('System noise (SD)');
zlabel('Coherence difference');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
title('coh(x,sysSUM(x,y)) > coh(x,sysID(x))');
subplot(1,2,2);
contour(X,Y,Z);
xlabel('Stimulus noise (SD)');
ylabel('System noise (SD)');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
c = colorbar;
c.Label.String = 'Coherence difference';
colormap(crameri('vik'));

temp = cellfun(@minus,coh_all_tbl.x_sSumxy_py,coh_all_tbl.x_sIdx_py,...
               'UniformOutput',false);
Z = reshape(cellfun(@mean,temp),numel(noise_stim),numel(noise_syst));
E = reshape(cellfun(@std,temp),numel(noise_stim),numel(noise_syst));
figure();
subplot(1,2,1);
surf(X,Y,Z); hold on
plot3([X(:),X(:)]', [Y(:),Y(:)]', [-E(:),E(:)]'+Z(:)','Color','k',...
      'LineWidth',1.5);
hold off
legend('mean','SD')
xlabel('Stimulus noise (SD)');
ylabel('System noise (SD)');
zlabel('Coherence difference');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
title('coh(x,sysSUM(x,y)|part y) > coh(x,sysID(x)|part y)');
subplot(1,2,2);
contour(X,Y,Z);
xlabel('Stimulus noise (SD)');
ylabel('System noise (SD)');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
c = colorbar;
c.Label.String = 'Coherence difference';
colormap(crameri('vik'));