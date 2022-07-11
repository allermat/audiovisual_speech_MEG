%% Simulation 1: Coherence between a single sinusoid and phase shifted set of sinusoids
clearvars;

% Generate a signal consisting of a sinusoid signal
fs = 256;                    % Sampling frequency (samples per second)
dt = 1/fs;                   % seconds per sample
StopTime = 5;                % length of signal in seconds
t = (0:dt:StopTime-dt)';     % seconds
F = 4;                       % Frequency of sinusoid (Hz)
% Generate sinusoid signal with the required frequency
sx = sin(2*pi*F*t);
% Generate 100 sinusoid signals with frequency F and randomly set phase
phy = 0.1*pi*randn(100,1);
sy = arrayfun(@(p) sin((2*pi*F*t)+p),phy,'UniformOutput',false);
phz = 0.3*pi*randn(100,1);
sz = arrayfun(@(p) sin((2*pi*F*t)+p),phz,'UniformOutput',false);
% Add random noise to the signals
x = sx;
y = cellfun(@(c) c + randn(fs*StopTime,1),sy,'UniformOutput',false);
y = cat(2,y{:});
z = cellfun(@(c) c + randn(fs*StopTime,1),sz,'UniformOutput',false);
z = cat(2,z{:});

% Estimate the magnitude-squared coherence of x and y
[pxy,fc] = cpsd(x,y,[],[],size(x,1),fs);
pxx = cpsd(x,x,[],[],size(x,1),fs);
pyy = cpsd(y,y,[],[],size(x,1),fs);
pxz = cpsd(x,z,[],[],size(x,1),fs);
pzz = cpsd(z,z,[],[],size(x,1),fs);

funCoh = @(xy,xx,yy) abs(mean(xy,2)).^2./(abs(mean(xx,2)).*abs(mean(yy,2)));
cxy = funCoh(pxy,pxx,pyy);
cxz = funCoh(pxz,pxx,pzz);

% Plot results
figure(); 
subplot(3,3,1:3);
plot(t,x(:,1)); hold on
plot(t,y(:,1));
plot(t,z(:,1));
xlabel('time (s)');
ylabel('signal ampl [a.u.]');
legend('signal x','signal y1','signal z1');
title('4 Hz sinusoids with added Gaussian noise')
h = get(gca,'Children');
set(gca,'Children',[h(3) h(2) h(1)])

subplot(3,3,4);
r = ones(size(phy,1),1);
[u,v] = pol2cart(phy,r);
compass(u,v);
title('Phase angles (y)');

subplot(3,3,5:6);
plot(fc,cxy);
xlim([0,40]);
ylim([0,1]);
xlabel('frequency (Hz)');
ylabel('coherence');
title('Coherence between x and y');

subplot(3,3,7);
r = ones(size(phz,1),1);
[u,v] = pol2cart(phz,r);
compass(u,v);
title('Phase angles (z)');

subplot(3,3,8:9);
plot(fc,cxz);
xlim([0,40]);
ylim([0,1]);
xlabel('frequency (Hz)');
ylabel('coherence');
title('Coherence between x and z');

%% Simulation 2: Generating two sets of signals coherent with each-other
% Generate a signal consisting of a sinusoid signal. 
fs = 256;                    % Sampling frequency (samples per second)
dt = 1/fs;                   % seconds per sample
StopTime = 5;                % length of signal in seconds
t = (0:dt:StopTime-dt)';     % seconds
F = 4;                       % Frequency of sinusoid (Hz)
% Generate 100 sinusoid signals with base frequency F and randomly set 
% amplitude and phase
ax = rand(100,1);
phx = 0.15*pi*randn(100,1);
sx = arrayfun(@(a,p) a*sin((2*pi*F*t)+p),ax,phx,'UniformOutput',false);
% Generate 100 sinusoid signals with base frequency F and randomly set 
% amplitude and phase
ay = rand(100,1);
phy = 0.15*pi*randn(100,1);
sy = arrayfun(@(a,p) a*sin((2*pi*F*t)+p),ay,phy,'UniformOutput',false);
% Add Gaussian noise to the sinusoids
x = cellfun(@(c) c + randn(fs*StopTime,1),sx,'UniformOutput',false);
x = cat(2,x{:});
y = cellfun(@(c) c + randn(fs*StopTime,1),sy,'UniformOutput',false);
y = cat(2,y{:});
% x = cat(2,sx{:});
% y = cat(2,sy{:});

% Estimate the magnitude-squared coherence of x and y
[cxy,fc] = coherence(x,y,fs);

% Plot results
figure(); 
subplot(3,4,1:4);
plot(t,x(:,1)); hold on
plot(t,y(:,1));
xlabel('time (s)');
ylabel('signal ampl [a.u.]');
legend('signal x1', 'signal y1');
title('Example 4 Hz signals with Gaussian noise added');

subplot(3,4,6);
[u,v] = pol2cart(phx,ax);
compass(u,v);
title('Phase angles (x)');
subplot(3,4,7);
[u,v] = pol2cart(phy,ay);
compass(u,v);
title('Phase angles (y)');

subplot(3,4,9:12);
plot(fc,cxy);
xlim([0,40]);
ylim([0,1]);
xlabel('frequency (Hz)');
ylabel('coherence');
title('Coherence between x and y');

%% Simulation 3: a system responding to coherent signals
% System respnse function: Output is the input plus noise
funSyst = @(x) x + 4*randn(size(x));
% System response to signals x and y separately
srx = funSyst(x);
sry = funSyst(y);
% Coherence between signal x and system response to x
[cxsrx,fc] = coherence(x,srx,fs);
% Coherence between signal y and system response to y
cysry = coherence(y,sry,fs);
% Coherence between the system's response and the other signal
cxsry = coherence(x,sry,fs);
cysrx = coherence(y,srx,fs);

% Partialize y out
cxsrx_py = coherence(x,srx,fs,'partialize',y);
% Partialize x out
cysry_px = coherence(y,sry,fs,'partialize',x);

% Coherence between the system's response and the other signal,
% partializing the signal to which the system responds
cxsry_py = coherence(x,sry,fs,'partialize',y);
% Partialize x out
cysrx_px = coherence(y,srx,fs,'partialize',x);

% Partialize a random signal out
cxsrx_prand = coherence(x,srx,fs,'partialize',0.01*randn(size(x)));
cysry_prand = coherence(y,sry,fs,'partialize',0.01*randn(size(y)));

figure(); 
subplot(4,4,1:4);
plot(t,x(:,1)); hold on
plot(t,y(:,1));
xlabel('time (s)');
ylabel('signal ampl [a.u.]');
legend('signal x1', 'signal y1');
title('Example 4 Hz signals with Gaussian noise added');

subplot(4,4,5:8);
plot(t,srx(:,1)); hold on
plot(t,sry(:,1));
xlabel('time (s)');
ylabel('signal ampl [a.u.]');
legend('response to x1', 'response to y1');
title('Example system responses to x and y');

subplot(4,4,9:10);
plot(fc,cxsrx); hold on;
plot(fc,cxsrx_py);
plot(fc,cysry);
plot(fc,cysry_px);
plot(fc,cxsry);
plot(fc,cysrx);
plot(fc,cxsry_py);
plot(fc,cysrx_px);
xlim([F-2,F+6]);
ylim([0,1]);
xlabel('frequency (Hz)');
ylabel('coherence');
legendStr = {'coh(x,sys(x))','coh(x,sys(x)|py)',...
             'coh(y,sys(y))','coh(y,sys(y)|px)',...
             'coh(x,sys(y))','coh(y,sys(x))',...
             'coh(x,sys(y)|py)','coh(y,sys(x)|px)'};
legend(legendStr);
title(sprintf('Coherence between signals and the system''s response'));

subplot(4,4,11:12);
hzIdx = fc == F; 
bar([cxsrx(hzIdx),cxsrx_py(hzIdx),cxsrx_prand(hzIdx),cxsry(hzIdx),cxsry_py(hzIdx);...
     cysry(hzIdx),cysry_px(hzIdx),cysry_prand(hzIdx),cysrx(hzIdx),cysrx_px(hzIdx)]);
ylabel('coherence');
set(gca,'XTickLabel',{'signal x','signal y'});
xtickangle(45);
legend({'sys(same)','sys(same)|p(other)','sys(same)|p(rand)','sys(other)',...
        'sys(other)|p(other)'},'Location','southeast');
title(sprintf('Coherence at %d Hz',F));

% Coherence with a system responding to both signals x and y
% Define system respnse function to two inputs: Output is the sum of the
% two inputs plus noise
funSystDual = @(x,y) x + y + 4*randn(size(x));
% funSystDual = @(x,y) (x .* y) + 4*randn(size(x));
% Compute system response to signals x and y
srd = funSystDual(x,y);

% Coherence between signal x and system response to x
[cxsrd,fc] = coherence(x,srd,fs);
% Coherence between signal y and system response to y
cysrd = coherence(y,srd,fs);

% Partialize y out
cxsrd_py = coherence(x,srx,fs,'partialize',y);
% Partialize x out
cysrd_px = coherence(y,sry,fs,'partialize',x);

% Continuing the previous figure
subplot(4,4,13:14);
plot(fc,cxsrd); hold on;
plot(fc,cxsrd_py);
plot(fc,cysrd);
plot(fc,cysrd_px);
xlim([F-2,F+6]);
ylim([0,1]);
xlabel('frequency (Hz)');
ylabel('coherence');
legendStr = {'coh(x,syst(x,y))','coh(x,syst(x,y)|party)',...
             'coh(y,syst(x,y))','coh(y,syst(x,y)|partx)'};
legend(legendStr);
title(sprintf('Coherence between signals and the system''s response'));

subplot(4,4,15:16);
hzIdx = fc == F; 
bar([cxsrd(hzIdx),cxsrd_py(hzIdx);...
     cysrd(hzIdx),cysrd_px(hzIdx)]);
ylabel('coherence');
set(gca,'XTickLabel',{'signal x','signal y'});
xtickangle(45);
legend({'coh','pcoh'},'Location','southeast');
title(sprintf('Coherence at %d Hz',F));

%% Simulation 4: Performing the above 100 times
nReps = 100;
% System response to signals x and y
srx_c = arrayfun(@(i) funSyst(x),1:nReps,'UniformOutput',false);
sry_c = arrayfun(@(i) funSyst(y),1:nReps,'UniformOutput',false);

% Coherence between signal x and system response to x
cxsrx_c = cellfun(@(i) coherence(x,i,fs),srx_c,'UniformOutput',false);
% Coherence between signal y and system response to y
cysry_c = cellfun(@(i) coherence(y,i,fs),sry_c,'UniformOutput',false);
% Partialize y out
cxsrx_py_c = cellfun(@(i) coherence(x,i,fs,'partialize',y),srx_c,'UniformOutput',false);
% Partialize x out
cysry_px_c = cellfun(@(i) coherence(y,i,fs,'partialize',x),sry_c,'UniformOutput',false);
% Partialize a random signal out
cxsrx_prand_c = cellfun(@(i) coherence(x,i,fs,'partialize',0.01*randn(size(x))),...
                        srx_c,'UniformOutput',false);
cysry_prand_c = cellfun(@(i) coherence(y,i,fs,'partialize',0.01*randn(size(y))),...
                        sry_c,'UniformOutput',false);

% Coherence between signal x and system response to y
cxsry_c = cellfun(@(i) coherence(x,i,fs),sry_c,'UniformOutput',false);
% Coherence between signal y and system response to x
cysrx_c = cellfun(@(i) coherence(y,i,fs),srx_c,'UniformOutput',false);
% Partialize y out
cxsry_py_c = cellfun(@(i) coherence(x,i,fs,'partialize',y),sry_c,'UniformOutput',false);
% Partialize x out
cysrx_px_c = cellfun(@(i) coherence(y,i,fs,'partialize',x),srx_c,'UniformOutput',false);

% System response to signals x and y
srd_c = arrayfun(@(i) funSystDual(x,y),1:nReps,'UniformOutput',false);

% Coherence between signal x and system response to x
cxsrd_c = cellfun(@(i) coherence(x,i,fs),srd_c,'UniformOutput',false);
% Coherence between signal y and system response to y
cysrd_c = cellfun(@(i) coherence(y,i,fs),srd_c,'UniformOutput',false);

% Partialize y out
cxsrd_py_c = cellfun(@(i) coherence(x,i,fs,'partialize',y),srd_c,'UniformOutput',false);
% Partialize x out
cysrd_px_c = cellfun(@(i) coherence(y,i,fs,'partialize',x),srd_c,'UniformOutput',false);

% computing mean and standard deviation across samples
cxsrx_m = mean(cat(2,cxsrx_c{:}),2);
cxsrx_std = std(cat(2,cxsrx_c{:}),[],2);
cysry_m = mean(cat(2,cysry_c{:}),2);
cysry_std = std(cat(2,cysry_c{:}),[],2);
cxsrx_py_m = mean(cat(2,cxsrx_py_c{:}),2);
cxsrx_py_std = std(cat(2,cxsrx_py_c{:}),[],2);
cysry_px_m = mean(cat(2,cysry_px_c{:}),2);
cysry_px_std = std(cat(2,cysry_px_c{:}),[],2);
cxsrx_prand_m = mean(cat(2,cxsrx_prand_c{:}),2);
cxsrx_prand_std = std(cat(2,cxsrx_prand_c{:}),[],2);
cysry_prand_m = mean(cat(2,cysry_prand_c{:}),2);
cysry_prand_std = std(cat(2,cysry_prand_c{:}),[],2);

cxsry_m = mean(cat(2,cxsry_c{:}),2);
cxsry_std = std(cat(2,cxsry_c{:}),[],2);
cysrx_m = mean(cat(2,cysrx_c{:}),2);
cysrx_std = std(cat(2,cysrx_c{:}),[],2);
cxsry_py_m = mean(cat(2,cxsry_py_c{:}),2);
cxsry_py_std = std(cat(2,cxsry_py_c{:}),[],2);
cysrx_px_m = mean(cat(2,cysrx_px_c{:}),2);
cysrx_px_std = std(cat(2,cysrx_px_c{:}),[],2);

cxsrd_m = mean(cat(2,cxsrd_c{:}),2);
cxsrd_std = std(cat(2,cxsrd_c{:}),[],2);
cysrd_m = mean(cat(2,cysrd_c{:}),2);
cysrd_std = std(cat(2,cysrd_c{:}),[],2);
cxsrd_py_m = mean(cat(2,cxsrd_py_c{:}),2);
cxsrd_py_std = std(cat(2,cxsrd_py_c{:}),[],2);
cysrd_px_m = mean(cat(2,cysrd_px_c{:}),2);
cysrd_px_std = std(cat(2,cysrd_px_c{:}),[],2);

figure(); 
subplot(2,2,1);
errorbar(fc,cxsrx_m,cxsrx_std); hold on;
errorbar(fc,cxsrx_py_m,cxsrx_py_std);
errorbar(fc,cysry_m,cysry_std);
errorbar(fc,cysry_px_m,cysry_px_std);
errorbar(fc,cxsry_m,cxsry_std);
errorbar(fc,cxsry_py_m,cxsry_py_std);
errorbar(fc,cysrx_m,cysrx_std);
errorbar(fc,cysrx_px_m,cysrx_px_std);
xlim([F-2,F+6]);
ylim([0,1]);
xlabel('frequency (Hz)');
ylabel('coherence');
legendStr = {'coh(x,sys(x))','coh(x,sys(x)|py)',...
             'coh(y,sys(y))','coh(y,sys(y)|px)',...
             'coh(x,sys(y))','coh(x,sys(y)|py)',...
             'coh(y,sys(x))','coh(y,sys(x)|px)'};
legend(legendStr);
title(sprintf('Coherence between signals and the system''s response \nmean+/-SD'));

subplot(2,2,2);
hzIdx = fc == F; 
barwitherr([cxsrx_std(hzIdx),cxsrx_py_std(hzIdx),cxsrx_prand_std(hzIdx),...
            cxsry_std(hzIdx),cxsry_py_std(hzIdx);...
            cysry_std(hzIdx),cysry_px_std(hzIdx),cysry_prand_std(hzIdx),...
            cysrx_std(hzIdx),cysrx_px_std(hzIdx)],...
           [cxsrx_m(hzIdx),cxsrx_py_m(hzIdx),cxsrx_prand_m(hzIdx),...
            cxsry_m(hzIdx),cxsry_py_m(hzIdx);...
            cysry_m(hzIdx),cysry_px_m(hzIdx),cysry_prand_m(hzIdx),...
            cysrx_m(hzIdx),cysrx_px_m(hzIdx)]);
ylabel('coherence');
set(gca,'XTickLabel',{'signal x','signal y'});
xtickangle(45);
legend({'sys(same)','sys(same)|p(other)','sys(same)|p(rand)','sys(other)',...
        'sys(other)|p(other)'},'Location','southeast');
title(sprintf('Coherence at %d Hz, mean+/-SD',F));

subplot(2,2,3);
errorbar(fc,cxsrd_m,cxsrd_std); hold on;
errorbar(fc,cxsrd_py_m,cxsrd_py_std);
errorbar(fc,cysrd_m,cysrd_std);
errorbar(fc,cysrd_px_m,cysrd_px_std);

xlim([F-2,F+6]);
ylim([0,1]);
xlabel('frequency (Hz)');
ylabel('coherence');
legendStr = {'coh(x,syst(x,y))','coh(x,syst(x,y)|party)',...
             'coh(y,syst(x,y))','coh(y,syst(x,y)|partx)'};
legend(legendStr);
title(sprintf('Coherence between signals and the system''s response \nmean+/-SD'));

subplot(2,2,4);
hzIdx = fc == F; 
barwitherr([cxsrd_std(hzIdx),cxsrd_py_std(hzIdx);...
            cysrd_std(hzIdx),cysrd_px_std(hzIdx)],...
           [cxsrd_m(hzIdx),cxsrd_py_m(hzIdx);...
            cysrd_m(hzIdx),cysrd_px_m(hzIdx)]);
ylabel('coherence');
set(gca,'XTickLabel',{'signal x','signal y'});
xtickangle(45);
legend({'coh','pcoh'},'Location','southeast');
title(sprintf('Coherence at %d Hz, mean+/-SD',F));
