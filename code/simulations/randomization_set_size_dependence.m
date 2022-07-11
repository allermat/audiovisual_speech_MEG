%% Simulation 1
% This is to show that if we randomly permute the elements of the vector,
% the number of elements assigned to their original position within the
% vector is on average the same no matter how long the vector is. However, 
% the proportion of these elements are much higher if the vector is
% shorter. This is somewhat common sense. The implication is that if we
% compute the coherence between the original vector and the randomly 
% permuted, we will observe higher coherence if the vector is shorter, 
% because there is higher proportion of elements stay at their original 
% place. 

clearvars;

nSim = 100;
nSubj = 14;
nRand = 100;
a = (1:10)';
b = (1:5)';
[sub_mean_num_identical_a, sub_mean_num_identical_b] = deal(zeros(nSim,nSubj));
[sub_mean_p_identical_a, sub_mean_p_identical_b] = deal(zeros(nSim,nSubj));
for iSim = 1:nSim
    for i = 1:nSubj
        p_identical_a = zeros(nRand,1);
        p_identical_b = zeros(nRand,1);
        num_identical_a = zeros(nRand,1);
        num_identical_b = zeros(nRand,1);
        for j = 1:nRand
            a_perm = a(randperm(numel(a)));
            b_perm = a(randperm(numel(b)));
            num_identical_a(j) = sum(a == a_perm);
            num_identical_b(j) = sum(b == b_perm);
            p_identical_a(j) = mean(a == a_perm);
            p_identical_b(j) = mean(b == b_perm);
        end
        sub_mean_num_identical_a(iSim,i) = mean(num_identical_a);
        sub_mean_num_identical_b(iSim,i) = mean(num_identical_b);
        sub_mean_p_identical_a(iSim,i) = mean(p_identical_a);
        sub_mean_p_identical_b(iSim,i) = mean(p_identical_b);
    end
end

figure();
subplot(2,2,1);
boxplot([mean(sub_mean_num_identical_a,2),mean(sub_mean_num_identical_b,2)]);
set(gca,'XTickLabel',{sprintf('N = %d', max(a)), sprintf('N = %d', max(b))})
ylabel('Number of identical elements');
title(sprintf('Average number of identical elements\nover 100 randomizations'));
subplot(2,2,2);
boxplot([std(sub_mean_num_identical_a,[],2),std(sub_mean_num_identical_b,[],2)]);
set(gca,'XTickLabel',{sprintf('N = %d', max(a)), sprintf('N = %d', max(b))})
ylabel('STD');
title(sprintf('STD of number of identical elements\nover 100 randomizations'));
subplot(2,2,3);
boxplot([mean(sub_mean_p_identical_a,2),mean(sub_mean_p_identical_b,2)]);
set(gca,'XTickLabel',{sprintf('N = %d', max(a)), sprintf('N = %d', max(b))})
ylabel('Proportion identical elements');
title(sprintf('Average proportion of identical elements\nover 100 randomizations'));
subplot(2,2,4);
boxplot([std(sub_mean_p_identical_a,[],2),std(sub_mean_p_identical_b,[],2)]);
set(gca,'XTickLabel',{sprintf('N = %d', max(a)), sprintf('N = %d', max(b))})
ylabel('STD');
title(sprintf('STD of proportion of identical elements\nover 100 randomizations'));

%% Simulation 2
% Single subject simulation, 100 samples vs all possible permutations. 
% This is to show that no matter how many permutations we perform, the 
% results of Simulation 1 stand. Specifically we perfomr all possible
% permutations and show that this is still the case. 
% The implication is that increasing the number of permutations to generate
% the permuted coherence in our analyses will not solve this issue. 
nSample = 100;
a = 10;
b = 5;
a_all = perms(1:a);
a_sampl = a_all(randsample(factorial(a), nSample), :);
b_all = perms(1:b);
b_sampl = b_all(randsample(factorial(b), nSample), :);

figure();
subplot(1,2,1);
boxplot(cat(1,mean(a_sampl == repmat(1:a,nSample,1),2), ...
            mean(b_sampl == repmat(1:b,nSample,1),2)), ...
        cat(1,ones(nSample,1),2*ones(nSample,1)));
set(gca,'XTickLabel',{sprintf('N = %d', a), sprintf('N = %d', b)})
ylabel('Porportion of identical elements');
title(sprintf('Average proportion of identical elements\nover %d randomizations',...
    nSample));
subplot(1,2,2);
boxplot(cat(1,mean(a_all == repmat(1:a,factorial(a),1),2), ...
            mean(b_all == repmat(1:b,factorial(b),1),2)), ...
        cat(1,ones(factorial(a),1),2*ones(factorial(b),1)));
set(gca,'XTickLabel',{sprintf('N = %d', a), sprintf('N = %d', b)})
ylabel('Porportion of identical elements');
title(sprintf('Average proportion of identical elements\nover all possible randomizations'));
mean(mean(a_all == repmat(1:a,factorial(a),1),2))
mean(mean(b_all == repmat(1:b,factorial(b),1),2))
% figure();
% subplot(1,2,1);
% histogram(mean(a_sampl == repmat(1:a,nSample,1),2));
% subplot(1,2,2);
% histogram(mean(b_sampl == repmat(1:b,nSample,1),2));

%% Simulation 3
% One possible solution to the above problem is if we only take those
% permutations where neither of the elements remain in their original
% position. This is called a dreangement, see: 
% https://en.wikipedia.org/wiki/Derangement. Surprisingly, as we will see,
% even if we use derangements to compute the permuted coherence, the
% dependence on the number of trials remains (i.e. smaller number of trials
% results higher permuted coherence). 
% The implication of these results is that it is imperative that we only
% compare coherence values RELATIVE to permuted across conditions,
% otherwise the results might be biased by the number of trials (if the
% number of trials across conditions is not uniform. Furthermore, it seems
% that with derangement the permuted coherence is slightly smaller
% compared to the ordinary permutation (which is expected). I will consider
% using the deranged coherence in place of the permuted coherence. I was
% hesitent doing this before as I was afraid that this is a biased sample
% of the original data, and indeed if I was performing a permutation test
% on the single subject level comparing true and permuted coherence, this 
% would bias the null-distribution, see Nichols and Holmes 2003, "Single
% Voxel Example" section. In our study however we do the statistical test 
% on the second, between subjects level, so generating our null condition
% using the deranged data seems a resonable choice. 
fs = 256;                    % Sampling frequency (samples per second)
StopTime = 5;                % length of signal in seconds
x = randn(fs*StopTime,100);
y = randn(fs*StopTime,50);

% Estimate the magnitude-squared coherence of x and y
nRand = 100;
[cxx,fc] = coherence(x,x,fs);
cx_randx = arrayfun(@(i) coherence(x,x(:,randperm(size(x,2))),fs),1:nRand,...
    'UniformOutput',false);
cy_randy = arrayfun(@(i) coherence(y,y(:,randperm(size(y,2))),fs),1:nRand,...
    'UniformOutput',false);
cx_randfullx = arrayfun(@(i) coherence(x,x(:,randpermfull(size(x,2))),fs),1:nRand,...
    'UniformOutput',false);
cy_randfully = arrayfun(@(i) coherence(y,y(:,randpermfull(size(y,2))),fs),1:nRand,...
    'UniformOutput',false);

% Organise data in a table
var_freq = repmat(fc,4,1);
var_set_size = repmat([100,50,100,50],size(fc,1),1);
var_set_size = categorical(var_set_size(:));
var_perm = repmat({'permute','permute','derange','derange'},size(fc,1),1);
var_perm = categorical(var_perm(:));
data = table(var_freq,var_set_size,var_perm,...
             cat(1,mean(cat(2,cx_randx{:}),2),mean(cat(2,cy_randy{:}),2),...
                 mean(cat(2,cx_randfullx{:}),2),mean(cat(2,cy_randfully{:}),2)),...
             'VariableNames',{'freq','set_size','perm_method','coh'});

% Plot results
% Plot across frequencies
figure();
plot(fc,data.coh(data.set_size == '100' & data.perm_method == 'permute')); hold on;
plot(fc,data.coh(data.set_size == '50' & data.perm_method == 'permute'));
plot(fc,data.coh(data.set_size == '100' & data.perm_method == 'derange'));
plot(fc,data.coh(data.set_size == '50' & data.perm_method == 'derange'));
legend({'N = 100, permutation', 'N = 50, permutation', ...
        'N = 100, derangement', 'N = 50, derangement'});
xlabel('frequency (Hz)');
ylabel('coherence');
title('Coherence between signal and permuted signal');
% Plot collapsed across frequencies
g = gramm('x',data.set_size,'y',data.coh,'color',data.perm_method);
g.stat_boxplot();
g.set_title('Coherence collapsed across frequencies');
g.set_names('x','set size','y','coherence','color','');
figure();
g.draw(false);
%% Simulation 4 
% Following from Simulation 3 it turned out that the dependence of
% permuted/deranged coherence on the number of trials follows from much
% more basic mathematical considerations. When computing coherence, the
% cross-spectra are averaged across trials and then the squared length of 
% the resultant vector is taken which is scaled by the powers of the two
% signals. The first step proved to be crucial: the more vectors (with 
% angles pointing randomly in any direction) we average, the shorter the
% resultant vector becomes. My friend, Dr. Gyula Toth has kindly provided a
% formal proof of that in the special case when the vectors are all unit
% length and the distribution of the angles is uniform. In this case the
% expected length of the resultant vector is 1/sqrt(N), where N is the
% number of vectors which are averaged. This simulation shows that this
% observation hold even if the vectors are not of unit length, but also
% random. Furthermore it is easy to change the code to show that the same
% is observed if the angles are normally distributed. 

% Perform a number of simulations
nRand = 500;
[x,y] = deal(cell(nRand,1));
% Set sizes
nx = 50;
ny = 100;
for j = 1:nRand
    % Generate unit power and random phase for 100 sinusoids
    ax = ones(nx,1);
    phx = 2*pi*rand(nx,1);
    x{j} = ax.*exp(1i.*phx);
    % Generate unit power and random phase for 50 sinusoids
    ay = ones(ny,1);
    phy = 2*pi*rand(ny,1);
    y{j} = ay.*exp(1i.*phy);
end
x_lengths = cellfun(@(c) abs(mean(c)),x);
x_angles = cellfun(@(c) angle(mean(c)),x);
y_lengths = cellfun(@(c) abs(mean(c)),y);
y_angles = cellfun(@(c) angle(mean(c)),y);

% Plot results
figure(); 
subplot(2,4,1);
[u,v] = pol2cart(angle(x{1}),abs(x{1}));
compass(u,v);
title(sprintf('Amplitudes and \nphase angles\n(N=%d)',nx));
subplot(2,4,2);
[u,v] = pol2cart(angle(y{1}),abs(y{1}));
compass(u,v,'r');
title(sprintf('Amplitudes and \nphase angles\n(N=%d)',ny));
subplot(2,4,5:6);
[u,v] = pol2cart(x_angles,x_lengths);
compass(u,v); hold on;
[u,v] = pol2cart(y_angles,y_lengths);
compass(u,v,'r'); 
title(sprintf('Mean amplitudes and phase angles\nacross %d simulations',nRand));
subplot(2,4,[3,4,7,8]);
boxplot(cat(2,x_lengths,y_lengths));
xlabel('Set size');
ylabel('Mean amplitude')
set(gca,'XTickLabel',{num2str(nx),num2str(ny)});
title(sprintf('Mean amplitudes across %d simulations',nRand));

%% Simulation 5
% Using the results of the previous simulations I generate two sets of
% sinusoids signals (similar to the ones in coherence_demo.m), one 
% including 50 and one including 100 signals. The question here is whether
% we can still see the dependenc of permuted and deranged coherence on
% sample size. Furthermore and most importantly, whether the difference 
% between true and permuted (resp. deranged) coherence is independent of
% sample size. I perform this simulation 100x. As the results show, both
% true-permuted and true-deranged show a tiny dependence on sample size.
% Hence we should match numbers of trials between compared conditions 
% (e.g., by subsampling)

% Sinusoid signal parameters
fs = 256;                    % Sampling frequency (samples per second)
dt = 1/fs;                   % seconds per sample
StopTime = 5;                % length of signal in seconds
t = (0:dt:StopTime-dt)';     % seconds
F = 4;                       % Frequency of sinusoid (Hz)
% System respnse function: Output is the input plus noise
funSyst = @(x) x + 4*randn(size(x));

% Generate 50 sinusoid signals with base frequency F and randomly set 
% amplitude and phase
nx = 50;
ax = rand(nx,1);
phx = 0.15*pi*randn(nx,1);
sx = arrayfun(@(a,p) a*sin((2*pi*F*t)+p),ax,phx,'UniformOutput',false);
% Generate 100 sinusoid signals with base frequency F and randomly set 
% amplitude and phase
ny = 100;
ay = rand(ny,1);
phy = 0.15*pi*randn(ny,1);
sy = arrayfun(@(a,p) a*sin((2*pi*F*t)+p),ay,phy,'UniformOutput',false);
% Add Gaussian noise to the sinusoids
x = cellfun(@(c) c + randn(fs*StopTime,1),sx,'UniformOutput',false);
x = cat(2,x{:});
y = cellfun(@(c) c + randn(fs*StopTime,1),sy,'UniformOutput',false);
y = cat(2,y{:});
% System response to signals x and y separately
srx = funSyst(x);
sry = funSyst(y);
% Estimate the magnitude-squared coherence between signals and system
% responses
% True coherence
[cx_srx,fc] = coherence(x,srx,fs);
cy_sry = coherence(y,sry,fs);
% Permuted/deranged coherence
nRand = 100;
cx_srx_r = arrayfun(@(i) coherence(x,funSyst(x(:,randperm(nx))),fs),1:nRand,...
                    'UniformOutput',false);
cx_srx_r = mean(cat(2,cx_srx_r{:}),2);
cy_sry_r = arrayfun(@(i) coherence(y,funSyst(y(:,randperm(ny))),fs),1:nRand,...
                    'UniformOutput',false);
cy_sry_r = mean(cat(2,cy_sry_r{:}),2);
cx_srx_d = arrayfun(@(i) coherence(x,funSyst(x(:,randpermfull(nx))),fs),1:nRand,...
                    'UniformOutput',false);
cx_srx_d = mean(cat(2,cx_srx_d{:}),2);
cy_sry_d = arrayfun(@(i) coherence(y,funSyst(y(:,randpermfull(ny))),fs),1:nRand,...
                    'UniformOutput',false);
cy_sry_d = mean(cat(2,cy_sry_d{:}),2);

% Plot results
figure(); 
subplot(4,4,1:4);
plot(t,x(:,1)); hold on
plot(t,y(:,1));
xlabel('time (s)');
ylabel('signal ampl [a.u.]');
legend('signal x1', 'signal y1');
title('Example 4 Hz signals with Gaussian noise added');

subplot(4,4,6);
[u,v] = pol2cart(phx,ax);
compass(u,v);
title(sprintf('Phase angles\n(x, N=%d)',nx));
subplot(4,4,7);
[u,v] = pol2cart(phy,ay);
compass(u,v);
title(sprintf('Phase angles\n(y, N=%d)',ny));

subplot(4,4,9:10);
plot(fc,cx_srx); hold on;
plot(fc,cy_sry);
plot(fc,cx_srx_r);
plot(fc,cy_sry_r);
plot(fc,cx_srx_d);
plot(fc,cy_sry_d);
xlim([F-2,F+4]);
ylim([0,1]);
xlabel('frequency (Hz)');
ylabel('coherence');
legendStr = {sprintf('coh true (N=%d)',nx),sprintf('coh true (N=%d)',ny),...
             sprintf('coh perm (N=%d)',nx),sprintf('coh perm (N=%d)',ny),...
             sprintf('coh derange (N=%d)',nx),sprintf('coh derange (N=%d)',ny)};
legend(legendStr);
title(sprintf('Coherence between signals and the system''s response'));

subplot(4,4,11:12);
hzIdx = fc == F; 
bar([cx_srx(hzIdx),cy_sry(hzIdx);...
     cx_srx_r(hzIdx),cy_sry_r(hzIdx);...
     cx_srx_d(hzIdx),cy_sry_d(hzIdx)]);
ylabel('coherence');
set(gca,'XTickLabel',{'true','permute','derange'});
xtickangle(45);
legend(sprintf('N = %d',nx),sprintf('N = %d',ny))
title(sprintf('Coherence at %d Hz',F));

% Repeat the above 100x
nRep = 100;
[cx_srx,cy_sry,cx_srx_r,cy_sry_r,cx_srx_d,cy_sry_d] = deal(cell(nRep,1));
parfor iRep = 1:nRep
    funSyst = @(x) x + 4*randn(size(x));
    % Signal x, N=50
    ax = rand(nx,1);
    phx = 0.15*pi*randn(nx,1);
    sx = arrayfun(@(a,p) a*sin((2*pi*F*t)+p),ax,phx,'UniformOutput',false);
    % Signal y, N=100
    ay = rand(ny,1);
    phy = 0.15*pi*randn(ny,1);
    sy = arrayfun(@(a,p) a*sin((2*pi*F*t)+p),ay,phy,'UniformOutput',false);
    % Add Gaussian noise to the sinusoids
    x = cellfun(@(c) c + randn(fs*StopTime,1),sx,'UniformOutput',false);
    x = cat(2,x{:});
    y = cellfun(@(c) c + randn(fs*StopTime,1),sy,'UniformOutput',false);
    y = cat(2,y{:});
    % System response to signals x and y separately
    srx = funSyst(x);
    sry = funSyst(y);
    % True coherence
    [cx_srx{iRep},fc] = coherence(x,srx,fs);
    cy_sry{iRep} = coherence(y,sry,fs);
    % Permuted/deranged coherence
    temp = arrayfun(@(i) coherence(x,funSyst(x(:,randperm(nx))),fs),1:nRand,...
        'UniformOutput',false);
    cx_srx_r{iRep} = mean(cat(2,temp{:}),2);
    temp = arrayfun(@(i) coherence(y,funSyst(y(:,randperm(ny))),fs),1:nRand,...
        'UniformOutput',false);
    cy_sry_r{iRep} = mean(cat(2,temp{:}),2);
    temp = arrayfun(@(i) coherence(x,funSyst(x(:,randpermfull(nx))),fs),1:nRand,...
        'UniformOutput',false);
    cx_srx_d{iRep} = mean(cat(2,temp{:}),2);
    temp = arrayfun(@(i) coherence(y,funSyst(y(:,randpermfull(ny))),fs),1:nRand,...
        'UniformOutput',false);
    cy_sry_d{iRep} = mean(cat(2,temp{:}),2);
end
% computing mean and standard deviation across samples
m_cx_srx = mean(cat(2,cx_srx{:}),2);
sd_cx_srx = std(cat(2,cx_srx{:}),[],2);
m_cy_sry = mean(cat(2,cy_sry{:}),2);
sd_cy_sry = std(cat(2,cy_sry{:}),[],2);
m_cx_srx_r = mean(cat(2,cx_srx_r{:}),2);
sd_cx_srx_r = std(cat(2,cx_srx_r{:}),[],2);
m_cy_sry_r = mean(cat(2,cy_sry_r{:}),2);
sd_cy_sry_r = std(cat(2,cy_sry_r{:}),[],2);
m_cx_srx_d = mean(cat(2,cx_srx_d{:}),2);
sd_cx_srx_d = std(cat(2,cx_srx_d{:}),[],2);
m_cy_sry_d = mean(cat(2,cy_sry_d{:}),2);
sd_cy_sry_d = std(cat(2,cy_sry_d{:}),[],2);
temp = cellfun(@minus,cx_srx,cx_srx_r,'UniformOutput',false);
m_cx_tr_r = mean(cat(2,temp{:}),2);
sd_cx_tr_r = std(cat(2,temp{:}),[],2);
temp = cellfun(@minus,cx_srx,cx_srx_d,'UniformOutput',false);
m_cx_tr_d = mean(cat(2,temp{:}),2);
sd_cx_tr_d = std(cat(2,temp{:}),[],2);
temp = cellfun(@minus,cy_sry,cy_sry_r,'UniformOutput',false);
m_cy_tr_r = mean(cat(2,temp{:}),2);
sd_cy_tr_r = std(cat(2,temp{:}),[],2);
temp = cellfun(@minus,cy_sry,cy_sry_d,'UniformOutput',false);
m_cy_tr_d = mean(cat(2,temp{:}),2);
sd_cy_tr_d = std(cat(2,temp{:}),[],2);

subplot(4,4,13:14);
hzIdx = fc == F;
barwitherr([sd_cx_srx(hzIdx),sd_cy_sry(hzIdx);...
            sd_cx_srx_r(hzIdx),sd_cy_sry_r(hzIdx);...
            sd_cx_srx_d(hzIdx),sd_cy_sry_d(hzIdx)],...
           [m_cx_srx(hzIdx),m_cy_sry(hzIdx);...
            m_cx_srx_r(hzIdx),m_cy_sry_r(hzIdx);...
            m_cx_srx_d(hzIdx),m_cy_sry_d(hzIdx)]);
ylabel('coherence');
set(gca,'XTickLabel',{'true','permute','derange'});
xtickangle(45);
legend(sprintf('N = %d',nx),sprintf('N = %d',ny))
title(sprintf('Coherence at %d Hz, mean+/-SD',F));

subplot(4,4,15:16);
hzIdx = fc == F;
barwitherr([sd_cx_tr_r(hzIdx),sd_cy_tr_r(hzIdx);...
            sd_cx_tr_d(hzIdx),sd_cy_tr_d(hzIdx)],...
           [m_cx_tr_r(hzIdx),m_cy_tr_r(hzIdx);...
            m_cx_tr_d(hzIdx),m_cy_tr_d(hzIdx)]);
ylabel('coherence');
set(gca,'XTickLabel',{'true-perm','true-derange'});
xtickangle(45);
legend(sprintf('N = %d',nx),sprintf('N = %d',ny))
title(sprintf('Coherence at %d Hz, mean+/-SD',F));

% Plot results as boxplots using gramm
% Organise data in a table
var_coh = cat(1,cellfun(@(x) x(hzIdx),cx_srx),cellfun(@(x) x(hzIdx),cx_srx_r),...
                cellfun(@(x) x(hzIdx),cx_srx_d),cellfun(@(x) x(hzIdx),cy_sry),...
                cellfun(@(x) x(hzIdx),cy_sry_r),cellfun(@(x) x(hzIdx),cy_sry_d));
var_set_size = repmat([nx,nx,nx,ny,ny,ny],nRep,1);
var_set_size = categorical(var_set_size(:));
var_perm = repmat({'true','permute','derange','true','permute','derange'},...
                  nRep,1);
var_perm = categorical(var_perm(:));
data = table(var_set_size,var_perm,var_coh,'VariableNames',...
            {'set_size','perm','coh'});

var_coh = cat(1,cellfun(@(x) x(hzIdx),cx_srx) - cellfun(@(x) x(hzIdx),cx_srx_r),...
                cellfun(@(x) x(hzIdx),cx_srx) - cellfun(@(x) x(hzIdx),cx_srx_d),...
                cellfun(@(x) x(hzIdx),cy_sry) - cellfun(@(x) x(hzIdx),cy_sry_r),...
                cellfun(@(x) x(hzIdx),cy_sry) - cellfun(@(x) x(hzIdx),cy_sry_d));
var_set_size = repmat([nx,nx,ny,ny],nRep,1);
var_set_size = categorical(var_set_size(:));
var_perm = repmat({'true-permute','true-derange','true-permute','true-derange'},...
                  nRep,1);
var_perm = categorical(var_perm(:));
data_diff = table(var_set_size,var_perm,var_coh,'VariableNames',...
            {'set_size','perm','coh'});
        
% Plotting data with box plots
g(1,1) = gramm('x',data.perm,'y',data.coh,'color',data.set_size);
g(1,1).stat_boxplot();
g(1,2) = gramm('x',data_diff.perm,'y',data_diff.coh,'color',data_diff.set_size);
g(1,2).stat_boxplot();
g.set_title(sprintf('Coherence at %d Hz',F));
g.set_names('x','','y','coherence','color','set size');
figure();
g.draw();

saveallfigures()
