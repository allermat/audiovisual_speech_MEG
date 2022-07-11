% Simulations to determine stimulus and signal SNR to key coherence
% results
clearvars;
currPool = gcp('nocreate');
if isempty(currPool)
    parpool('local');
end
% Preparing the parameter space
% grid resolution
weight_x = 0:5;
weight_y = 0:5; 
% Parameter names
parameterNames = {'weight_x','weight_y'};
gridVectors = {weight_x,weight_y};
nParameters = numel(gridVectors);
% Full factorial expansion of the specified parameter vectors
coords = cell(1,nParameters);
[coords{:}] = ndgrid(gridVectors{:});
coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
% Matrix of all possible parameter combinations: nRows = number of
% combinations, columns correspond to parameters
paramCombinations = cat(2,coords{:});
% Remove the case where both stimulus weights are zero
% paramCombinations(1,:) = [];

% Generating two sets of signals coherent with each-other
% Generate a signal consisting of a sinusoid signal.
fs = 256;                    % Sampling frequency (samples per second)
dt = 1/fs;                   % seconds per sample
StopTime = 5;                % length of signal in seconds
t = (0:dt:StopTime-dt)';     % seconds
F = 4;                       % Frequency of sinusoid (Hz)
nSig = 100; 
nSubj = 100;

% Noise parameters are fixed
noise_stim = 1;
noise_syst = 4;

% Handy functions to compute mean and std across repetitions
fun_mean = @(x) mean(cat(2,x{:}),2);
fun_std = @(x) std(cat(2,x{:}),[],2);

% Array for storing outputs
coh_all = cell(size(paramCombinations,1),4);
parfor_progress(size(paramCombinations,1));
parfor iGrid = 2:size(paramCombinations,1)
    % Leaving the case out where both stimulus weights are zero
    
    wx = paramCombinations(iGrid,1); % Amount of noise added to the signals
    wy = paramCombinations(iGrid,2); % Amount of noise added by the system
       
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
    noise_x = noise_stim;
    noise_y = noise_x;
    x = cellfun(@(c) c + noise_x*randn(fs*StopTime,1),sx,'UniformOutput',false);
    x = cat(2,x{:});
    y = cellfun(@(c) c + noise_y*randn(fs*StopTime,1),sy,'UniformOutput',false);
    y = cat(2,y{:});
    
    % Define system responses 
    syst_id = @(x) x + noise_syst*randn(size(x));
    % I normalise by the sum of the weights to keep the amount of signal
    % relative to noise constant
    syst_sum = @(x,y,wx,wy) (wx*x + wy*y)/(wx+wy) + noise_syst*randn(size(x));
    
    % Compute system responses to input
    sr_sum = arrayfun(@(i) syst_sum(x,y,wx,wy),1:nSubj,'UniformOutput',false);
    % Estimate the magnitude-squared coherence
    [c_sumxy_x,fc] = cellfun(@(a) coherence(x,a,fs),sr_sum,'UniformOutput',false);
    fc = fc{1};
    c_sumxy_x_py = cellfun(@(a) coherence(x,a,fs,'partialize',y),sr_sum,...
                           'UniformOutput',false);
    c_sumxy_y = cellfun(@(a) coherence(y,a,fs),sr_sum,'UniformOutput',false);
    c_sumxy_y_px = cellfun(@(a) coherence(y,a,fs,'partialize',x),sr_sum,...
                           'UniformOutput',false);
    
    hzIdx = fc == F;
    temp = {cat(2,c_sumxy_x{:})',cat(2,c_sumxy_x_py{:})',...
            cat(2,c_sumxy_y{:})',cat(2,c_sumxy_y_px{:})'};
    coh_all(iGrid,:) = cellfun(@(x) x(:,hzIdx),temp,'UniformOutput',false);
    
    parfor_progress;
end
parfor_progress(0);
[coh_all{1,:}] = deal(NaN(nSubj,1));
% Plot figures
coh_all_tbl = cell2table(coh_all,'VariableNames',{'x_sSumxy','x_sSumxy_py',...
                         'y_sSumxy','y_sSumxy_px'});
X = reshape(paramCombinations(:,1),numel(weight_x),numel(weight_y));
Y = reshape(paramCombinations(:,2),numel(weight_x),numel(weight_y));
% Difference between partial coherences
temp = cellfun(@minus,coh_all_tbl.x_sSumxy_py,coh_all_tbl.y_sSumxy_px,...
               'UniformOutput',false);
Z = reshape(cellfun(@mean,temp),numel(weight_x),numel(weight_y));
E = reshape(cellfun(@std,temp),numel(weight_x),numel(weight_y));
figure();
subplot(1,2,1);
surf(X,Y,Z); hold on
plot3([X(:),X(:)]', [Y(:),Y(:)]', [-E(:),E(:)]'+Z(:)','Color','k',...
      'LineWidth',1.5);
hold off
legend('mean','SD')
xlabel('x weight');
ylabel('y weight');
zlabel('Coherence difference');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
title('coh(x,sysSUM(x,y)|part y) > coh(y,sysSUM(x,y)|part x)');
subplot(1,2,2);
contour(X,Y,Z); hold on
line([min(weight_x),max(weight_x)]',[min(weight_x),max(weight_x)]',...
     'Color',[0.5,0.5,0.5],'LineStyle','--');
xlabel('x weight');
ylabel('y weight');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
c = colorbar;
c.Label.String = 'Coherence difference';
colormap(crameri('vik'));
% Difference between coherences
temp = cellfun(@minus,coh_all_tbl.x_sSumxy,coh_all_tbl.y_sSumxy,...
               'UniformOutput',false);
Z = reshape(cellfun(@mean,temp),numel(weight_x),numel(weight_y));
E = reshape(cellfun(@std,temp),numel(weight_x),numel(weight_y));
figure();
subplot(1,2,1);
surf(X,Y,Z); hold on
plot3([X(:),X(:)]', [Y(:),Y(:)]', [-E(:),E(:)]'+Z(:)','Color','k',...
      'LineWidth',1.5);
hold off
legend('mean','SD')
xlabel('x weight');
ylabel('y weight');
zlabel('Coherence difference');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
title('coh(x,sysSUM(x,y)) > coh(y,sysSUM(x,y))');
subplot(1,2,2);
contour(X,Y,Z); hold on
line([min(weight_x),max(weight_x)]',[min(weight_x),max(weight_x)]',...
     'Color',[0.5,0.5,0.5],'LineStyle','--');
xlabel('x weight');
ylabel('y weight');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
c = colorbar;
c.Label.String = 'Coherence difference';
colormap(crameri('vik'));
% Partial coherence x
temp = coh_all_tbl.x_sSumxy_py;
Z = reshape(cellfun(@mean,temp),numel(weight_x),numel(weight_y));
E = reshape(cellfun(@std,temp),numel(weight_x),numel(weight_y));
figure();
subplot(1,2,1);
surf(X,Y,Z); hold on
plot3([X(:),X(:)]', [Y(:),Y(:)]', [-E(:),E(:)]'+Z(:)','Color','k',...
      'LineWidth',1.5);
hold off
legend('mean','SD')
xlabel('x weight');
ylabel('y weight');
zlabel('Coherence');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
title('coh(x,sysSUM(x,y)|part y)');
subplot(1,2,2);
contour(X,Y,Z); hold on
line([min(weight_x),max(weight_x)]',[min(weight_x),max(weight_x)]',...
     'Color',[0.5,0.5,0.5],'LineStyle','--');
xlabel('x weight');
ylabel('y weight');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
c = colorbar;
c.Label.String = 'Coherence';
colormap(crameri('vik'));
% Partial coherence y
temp = coh_all_tbl.y_sSumxy_px;
Z = reshape(cellfun(@mean,temp),numel(weight_x),numel(weight_y));
E = reshape(cellfun(@std,temp),numel(weight_x),numel(weight_y));
figure();
subplot(1,2,1);
surf(X,Y,Z); hold on
plot3([X(:),X(:)]', [Y(:),Y(:)]', [-E(:),E(:)]'+Z(:)','Color','k',...
      'LineWidth',1.5);
hold off
legend('mean','SD')
xlabel('x weight');
ylabel('y weight');
zlabel('Coherence');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
title('coh(y,sysSUM(x,y)|part x)');
subplot(1,2,2);
contour(X,Y,Z); hold on
line([min(weight_x),max(weight_x)]',[min(weight_x),max(weight_x)]',...
     'Color',[0.5,0.5,0.5],'LineStyle','--');
xlabel('x weight');
ylabel('y weight');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
c = colorbar;
c.Label.String = 'Coherence';
colormap(crameri('vik'));
% Coherence x
temp = coh_all_tbl.x_sSumxy;
Z = reshape(cellfun(@mean,temp),numel(weight_x),numel(weight_y));
E = reshape(cellfun(@std,temp),numel(weight_x),numel(weight_y));
figure();
subplot(1,2,1);
surf(X,Y,Z); hold on
plot3([X(:),X(:)]', [Y(:),Y(:)]', [-E(:),E(:)]'+Z(:)','Color','k',...
      'LineWidth',1.5);
hold off
legend('mean','SD')
xlabel('x weight');
ylabel('y weight');
zlabel('Coherence');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
title('coh(x,sysSUM(x,y)');
subplot(1,2,2);
contour(X,Y,Z); hold on
line([min(weight_x),max(weight_x)]',[min(weight_x),max(weight_x)]',...
     'Color',[0.5,0.5,0.5],'LineStyle','--');
xlabel('x weight');
ylabel('y weight');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
c = colorbar;
c.Label.String = 'Coherence';
colormap(crameri('vik'));
% Coherence y
temp = coh_all_tbl.y_sSumxy;
Z = reshape(cellfun(@mean,temp),numel(weight_x),numel(weight_y));
E = reshape(cellfun(@std,temp),numel(weight_x),numel(weight_y));
figure();
subplot(1,2,1);
surf(X,Y,Z); hold on
plot3([X(:),X(:)]', [Y(:),Y(:)]', [-E(:),E(:)]'+Z(:)','Color','k',...
      'LineWidth',1.5);
hold off
legend('mean','SD')
xlabel('x weight');
ylabel('y weight');
zlabel('Coherence');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
title('coh(y,sysSUM(x,y))');
subplot(1,2,2);
contour(X,Y,Z); hold on
line([min(weight_x),max(weight_x)]',[min(weight_x),max(weight_x)]',...
     'Color',[0.5,0.5,0.5],'LineStyle','--');
xlabel('x weight');
ylabel('y weight');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
c = colorbar;
c.Label.String = 'Coherence';
colormap(crameri('vik'));
% Difference betweeen coherence and partial coherence x
temp = cellfun(@minus,coh_all_tbl.x_sSumxy,coh_all_tbl.x_sSumxy_py,...
               'UniformOutput',false);
Z = reshape(cellfun(@mean,temp),numel(weight_x),numel(weight_y));
E = reshape(cellfun(@std,temp),numel(weight_x),numel(weight_y));
figure();
subplot(1,2,1);
surf(X,Y,Z); hold on
plot3([X(:),X(:)]', [Y(:),Y(:)]', [-E(:),E(:)]'+Z(:)','Color','k',...
      'LineWidth',1.5);
hold off
legend('mean','SD')
xlabel('x weight');
ylabel('y weight');
zlabel('Coherence difference');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
title('coh(x,sysSUM(x,y)) > coh(x,sysSUM(x,y)|part y)');
subplot(1,2,2);
contour(X,Y,Z); hold on
line([min(weight_x),max(weight_x)]',[min(weight_x),max(weight_x)]',...
     'Color',[0.5,0.5,0.5],'LineStyle','--');
xlabel('x weight');
ylabel('y weight');
caxis([-max(abs(Z(:))),max(abs(Z(:)))]);
c = colorbar;
c.Label.String = 'Coherence difference';
colormap(crameri('vik'));
