function [ftStat_out,ftDataGrAvg_out] = coh_sensor_group_stats(test,varargin)
% Sensor space coherence analysis statistical tests
%
% Copyright(C) Mate Aller 2020
% allermat@gmail.com

% Parsing input, checking matlab
p = inputParser;

validMegChannels = {'MEG','MEGGRAD','MEGMAG'};
validStats = {'tstat','zdiff'};
validTests = {'basicAudVis','aud_>_vis','accClar_x_visEnh'};

addRequired(p,'test',@(x) ismember(x,validTests));
addParameter(p,'nPerm',5000,@(x) validateattributes(x,{'double'},{'scalar',...
                                'integer','positive'}));
addParameter(p,'chanMeg','MEGMAG',@(x) ismember(x,validMegChannels));
addParameter(p,'combinePlanar',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'stat','tstat',@(x) ismember(x,validStats));
addParameter(p,'avgOverFreq',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,test,varargin{:});

test = p.Results.test;
nPerm = p.Results.nPerm;
chanMeg = p.Results.chanMeg;
combinePlanar = p.Results.combinePlanar;
stat = p.Results.stat;
avgOverFreq = p.Results.avgOverFreq;

% Selecting appropriate statfun
switch stat
    case 'zdiff', statfun = @statZdiff;
    case 'tstat', statfun = @statTtest;
end

% Define data folders
fprintf('Loading data...\n');
subjList = dir(fullfile(AVSM_setupdir('analysis_megcoherence'),'sub-*'));
subjList = {subjList.name}';
subjDirList = fullfile(AVSM_setupdir('analysis_megcoherence'),subjList);

switch test
    case 'basicAudVis'
        % File names to be loaded
        fileMatchStr = {'ftmeg_coh_aud-allAud';...
                        'ftmeg_coh_aud-allAud_perm';...
                        'ftmeg_coh_vis-allVis';...
                        'ftmeg_coh_vis-allVis_perm'};
        if strcmp(chanMeg,'MEGGRAD') && combinePlanar
            fileMatchStr = strcat(fileMatchStr,'_combPlanar');
        end
        % Loading data
        filesLoaded = cell(size(fileMatchStr));
        for i = 1:numel(fileMatchStr)
            temp = cellfun(@load,fullfile(subjDirList,...
                                          strcat(subjList,'_',fileMatchStr{i},'.mat')));
            filesLoaded{i} = {temp.ftData_coh}';
        end
        [ftDataCell_coh_aud,ftDataCell_coh_aud_perm,ftDataCell_coh_vis,...
            ftDataCell_coh_vis_perm] = deal(filesLoaded{:});
        clear filesLoaded temp
        
        % Statistical test settings
        cfg = [];
        cfg.channel = chanMeg;
        cfg.numrandomization = nPerm;
        cfg.frequency = [2,8];
        if avgOverFreq, cfg.avgoverfreq = 'yes'; else, cfg.avgoverfreq = 'no'; end
        % Computing stats
        ftStat_aud = statfun(ftDataCell_coh_aud,ftDataCell_coh_aud_perm,cfg);
        ftStat_vis = statfun(ftDataCell_coh_vis,ftDataCell_coh_vis_perm,cfg);

        ftStat_out = {ftStat_aud;ftStat_vis};
        
        if nargout > 2
            % Compute grand average topographies as well
            % Some additional tweaks to pretend the cohernece data are
            % powerspectra, so ft_freqgrandaverage doesn't crash
            ftDataCell_coh_aud = cellfun(@ftDataCohToFreq,ftDataCell_coh_aud,...
                'UniformOutput',false);
            ftDataCell_coh_vis = cellfun(@ftDataCohToFreq,ftDataCell_coh_vis,...
                'UniformOutput',false);
            cfg = struct();
            cfg.parameter = 'powspctrm';
            cfg.foilim = [2,8];
            ftDataGrAvg_coh_aud = ft_freqgrandaverage(cfg,ftDataCell_coh_aud{:});
            ftDataGrAvg_coh_vis = ft_freqgrandaverage(cfg,ftDataCell_coh_vis{:});
            ftDataGrAvg_out = {ftDataGrAvg_coh_aud;ftDataGrAvg_coh_vis};
        end
        
    case 'aud_>_vis'
        % Use only AV trials
        fileMatchStr = {'ftmeg_coh_aud-AV70';'ftmeg_coh_aud-AV20';...
                        'ftmeg_coh_vis-AV70';'ftmeg_coh_vis-AV20'};
        if strcmp(chanMeg,'MEGMAG') || (strcmp(chanMeg,'MEGGRAD') && combinePlanar)
            fileMatchStr = strcat(fileMatchStr,'_combPlanar');
        end
        % Loading data
        filesLoaded = cell(size(fileMatchStr));
        for i = 1:numel(fileMatchStr)
            temp = cellfun(@load,fullfile(subjDirList,strcat(fileMatchStr{i},'.mat')));
            filesLoaded{i} = {temp.ftData_coh}';
        end
        [ftDataCell_coh_aud_AV70,ftDataCell_coh_aud_AV20,...
         ftDataCell_coh_vis_AV70,ftDataCell_coh_vis_AV20] = deal(filesLoaded{:});
        clear filesLoaded temp
        
        % Statistical test settings
        cfg = [];
        cfg.channel = chanMeg;
        cfg.numrandomization = nPerm;
        if avgOverFreq, cfg.avgoverfreq = 'yes'; else, cfg.avgoverfreq = 'no'; end
        cfg.frequency = [2,6]; % Based on basic aud and vis coherence results
        % Computing stats
        ftDataCell_coh_aud = cellfun(@(x,y) ft_mathoperation(x,y,'mean'),...
            ftDataCell_coh_aud_AV20,ftDataCell_coh_aud_AV70,'UniformOutput',false);
        ftDataCell_coh_vis = cellfun(@(x,y) ft_mathoperation(x,y,'mean'),...
            ftDataCell_coh_vis_AV20,ftDataCell_coh_vis_AV70,'UniformOutput',false);
        cfg.tail = 1; % right-tailed test aud > vis
        ftStat_aud_vs_vis = statfun(ftDataCell_coh_aud,ftDataCell_coh_vis,cfg);
        cfg.tail = -1; % left-tailed test aud < vis
        ftStat_vis_vs_aud = statfun(ftDataCell_coh_aud,ftDataCell_coh_vis,cfg);
        ftStat_out = {ftStat_aud_vs_vis,ftStat_vis_vs_aud};
        
        if nargout > 2
            % Compute grand average topographies as well
            % Some additional tweaks to pretend the cohernece data are
            % powerspectra, so ft_freqgrandaverage doesn't crash
            ftDataCell_coh_aud = cellfun(@ftDataCohToFreq,ftDataCell_coh_aud,...
                'UniformOutput',false);
            ftDataCell_coh_vis = cellfun(@ftDataCohToFreq,ftDataCell_coh_vis,...
                'UniformOutput',false);
            cfg = struct();
            cfg.parameter = 'powspctrm';
            cfg.foilim = [2,9];
            ftDataGrAvg_coh_aud = ft_freqgrandaverage(cfg,ftDataCell_coh_aud{:});
            ftDataGrAvg_coh_vis = ft_freqgrandaverage(cfg,ftDataCell_coh_vis{:});
            ftDataGrAvg_out = {ftDataGrAvg_coh_aud;ftDataGrAvg_coh_vis};
        end
        
    % Main effects and interaction for accoustic clarity and visual
    % enhancement
    case 'accClar_x_visEnh'
        % File names to be loaded
        fileMatchStr = {'ftmeg_coh_aud-AO70';'ftmeg_coh_aud-AV70';...
                        'ftmeg_coh_aud-AO20';'ftmeg_coh_aud-AV20';...
                        'ftmeg_coh_vis-AO70';'ftmeg_coh_vis-AV70';...
                        'ftmeg_coh_vis-AO20';'ftmeg_coh_vis-AV20'};
        if strcmp(chanMeg,'MEGMAG') || (strcmp(chanMeg,'MEGGRAD') && combinePlanar)
            fileMatchStr = strcat(fileMatchStr,'_combPlanar');
        end
        % Loading data
        filesLoaded = cell(size(fileMatchStr));
        for i = 1:numel(fileMatchStr)
            temp = cellfun(@load,fullfile(subjDirList,strcat(fileMatchStr{i},'.mat')));
            filesLoaded{i} = {temp.ftData_coh}';
            % Converting coherency to coherence             
            % filesLoaded{i} = cellfun(@ftDataCoherency2coherence,filesLoaded{i},...
            %                          'UniformOutput',false);
        end
        [ftDataCell_coh_aud_AO70,ftDataCell_coh_aud_AV70,...
         ftDataCell_coh_aud_AO20,ftDataCell_coh_aud_AV20,...
         ftDataCell_coh_vis_AO70,ftDataCell_coh_vis_AV70,...
         ftDataCell_coh_vis_AO20,ftDataCell_coh_vis_AV20] = deal(filesLoaded{:});
        clear filesLoaded temp
        
        % Statistical test settings
        cfg = [];
        cfg.channel = chanMeg;
        cfg.numrandomization = nPerm;
        if avgOverFreq, cfg.avgoverfreq = 'yes'; else, cfg.avgoverfreq = 'no'; end
        
        % Main effect of accoustic clarity
        % Auditory coherence
        ftDataCell_coh_aud_70 = cellfun(@(x,y) ft_mathoperation(x,y,'mean'),...
            ftDataCell_coh_aud_AO70,ftDataCell_coh_aud_AV70,'UniformOutput',false);
        ftDataCell_coh_aud_20 = cellfun(@(x,y) ft_mathoperation(x,y,'mean'),...
            ftDataCell_coh_aud_AO20,ftDataCell_coh_aud_AV20,'UniformOutput',false);
        % Computing stats
        ftStat_aud_accClar = statfun(ftDataCell_coh_aud_70,ftDataCell_coh_aud_20,cfg);
        % Visual coherence
        ftDataCell_coh_vis_70 = cellfun(@(x,y) ft_mathoperation(x,y,'mean'),...
            ftDataCell_coh_vis_AO70,ftDataCell_coh_vis_AV70,'UniformOutput',false);
        ftDataCell_coh_vis_20 = cellfun(@(x,y) ft_mathoperation(x,y,'mean'),...
            ftDataCell_coh_vis_AO20,ftDataCell_coh_vis_AV20,'UniformOutput',false);
        % Computing stats
        ftStat_vis_accClar = statfun(ftDataCell_coh_vis_70,ftDataCell_coh_vis_20,cfg);
        
        % Main effect of visual enhancement
        % Auditory coherence
        ftDataCell_coh_aud_AV = cellfun(@(x,y) ft_mathoperation(x,y,'mean'),...
            ftDataCell_coh_aud_AV20,ftDataCell_coh_aud_AV70,'UniformOutput',false);
        ftDataCell_coh_aud_AO = cellfun(@(x,y) ft_mathoperation(x,y,'mean'),...
            ftDataCell_coh_aud_AO20,ftDataCell_coh_aud_AO70,'UniformOutput',false);
        % Computing stats
        ftStat_aud_visEnh = statfun(ftDataCell_coh_aud_AV,ftDataCell_coh_aud_AO,cfg);
        % Visual coherence
        ftDataCell_coh_vis_AV = cellfun(@(x,y) ft_mathoperation(x,y,'mean'),...
            ftDataCell_coh_vis_AV20,ftDataCell_coh_vis_AV70,'UniformOutput',false);
        ftDataCell_coh_vis_AO = cellfun(@(x,y) ft_mathoperation(x,y,'mean'),...
            ftDataCell_coh_vis_AO20,ftDataCell_coh_vis_AO70,'UniformOutput',false);
        % Computing stats
        ftStat_vis_visEnh = statfun(ftDataCell_coh_vis_AV,ftDataCell_coh_vis_AO,cfg);
                
        % Interaction
        % Auditory coherence
        ftDataCell_coh_aud_AVdiff70 = cellfun(@(x,y) ft_mathoperation(x,y,'subtract'),...
            ftDataCell_coh_aud_AV70,ftDataCell_coh_aud_AO70,'UniformOutput',false);
        ftDataCell_coh_aud_AVdiff20 = cellfun(@(x,y) ft_mathoperation(x,y,'subtract'),...
            ftDataCell_coh_aud_AV20,ftDataCell_coh_aud_AO20,'UniformOutput',false);
        % Computing stats
        ftStat_aud_int = statfun(ftDataCell_coh_aud_AVdiff70,ftDataCell_coh_aud_AVdiff20,cfg);
        % Visual coherence
        ftDataCell_coh_vis_AVdiff70 = cellfun(@(x,y) ft_mathoperation(x,y,'subtract'),...
            ftDataCell_coh_vis_AV70,ftDataCell_coh_vis_AO70,'UniformOutput',false);
        ftDataCell_coh_vis_AVdiff20 = cellfun(@(x,y) ft_mathoperation(x,y,'subtract'),...
            ftDataCell_coh_vis_AV20,ftDataCell_coh_vis_AO20,'UniformOutput',false);
        % Computing stats
        ftStat_vis_int = statfun(ftDataCell_coh_vis_AVdiff70,ftDataCell_coh_vis_AVdiff20,cfg);

        ftStat_out = {ftStat_aud_accClar;ftStat_vis_accClar;
                      ftStat_aud_visEnh;ftStat_vis_visEnh;...
                      ftStat_aud_int;ftStat_vis_int};
        
        if nargout > 2
            % Compute grand average topographies as well
            % Some additional tweaks to pretend the cohernece data are
            % powerspectra, so ft_freqgrandaverage doesn't crash
            ftDataCell_coh_aud = cellfun(@ftDataCohToFreq,ftDataCell_coh_aud,...
                'UniformOutput',false);
            ftDataCell_coh_vis = cellfun(@ftDataCohToFreq,ftDataCell_coh_vis,...
                'UniformOutput',false);
            cfg = struct();
            cfg.parameter = 'powspctrm';
            cfg.foilim = [2,9];
            ftDataGrAvg_coh_aud = ft_freqgrandaverage(cfg,ftDataCell_coh_aud{:});
            ftDataGrAvg_coh_vis = ft_freqgrandaverage(cfg,ftDataCell_coh_vis{:});
            ftDataGrAvg_out = {ftDataGrAvg_coh_aud;ftDataGrAvg_coh_vis};
        end
end

end


function ftData_out = ft_mathoperation(ftData_in1,ftData_in2,operation)
% Utility function to perform math operations on FT data
% 
% I know, there is ft_math for this purpose buth it is terribly slow

ftData_out = ftData_in1;
switch operation
    case 'mean'
        dim2cat = ndims(ftData_out.cohspctrm)+1;
        ftData_out.cohspctrm = mean(cat(dim2cat,...
            ftData_in1.cohspctrm,ftData_in2.cohspctrm),dim2cat);
    case 'subtract'
        ftData_out.cohspctrm = ftData_in1.cohspctrm - ftData_in2.cohspctrm;
    otherwise
        error('Unsopported operation');
end

end


function ftStat = statZdiff(ftDataCell_coh_c1,ftDataCell_coh_c2,cfg)

% I am using the Z-transformed coherence difference as the test statistic
% here with subjects as units of observation, for further info see: 
% https://mailman.science.ru.nl/pipermail/fieldtrip/2011-September/004279.html
ftDataCell_cohZ_c1_vs_c2 = cellfun(@(x,y) ftDataCohZtransform(x,y),ftDataCell_coh_c1,...
                                  ftDataCell_coh_c2,'UniformOutput',false);

% To compare the Z-transformed coherence differences to zero we need a set
% of null data
ftDataCell_coh_null = ftDataCell_coh_c1;
for i = 1:numel(ftDataCell_coh_null)
    ftDataCell_coh_null{i}.cohspctrm(:) = 0;
end
                          
% Some additional tweaks to pretend the cohernece data are powerspectra, so
% ft_freqstatistics doesn't crash, see: 
% https://mailman.science.ru.nl/pipermail/fieldtrip/2018-October/038311.html
ftDataCell_cohZ_c1_vs_c2 = cellfun(@ftDataCohToFreq,ftDataCell_cohZ_c1_vs_c2,...
                                  'UniformOutput',false);
ftDataCell_coh_null = cellfun(@ftDataCohToFreq,ftDataCell_coh_null,...
                              'UniformOutput',false);

% Set up frequency statistic parameters
cfg.channel = ft_getopt(cfg, 'channel', 'all');
combinePlanar = any(~cellfun(@isempty,regexp(ftDataCell_coh_c1{1}.label,'.*\+.*','once')));
if strcmp(cfg.channel,'MEGMAG')
    layout = 'neuromag306mag.lay';
elseif strcmp(cfg.channel,'MEGGRAD') && ~combinePlanar
    layout = 'neuromag306planar.lay';
elseif strcmp(cfg.channel,'MEGGRAD') && combinePlanar
    layout = 'neuromag306cmb.lay';
end
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesDiff'; 
cfg.numrandomization = ft_getopt(cfg,'numrandomization',5000);
cfg.correctm = 'cluster';
cfg.alpha = 0.05; 
cfg.tail = ft_getopt(cfg,'tail',1); 
if cfg.tail == 0,  cfg.correcttail = 'alpha'; end
% Use parametric thresholding as the statfun doesn't compute it and we
% anyway use a Z-score as a statistic
cfg.clustercritval = norminv(0.95); % Corresponding right tailed p = 0.05
cfg.frequency = ft_getopt(cfg,'frequency',[2,9]); % Based on auditory and lip envelope coherence
cfg.avgoverfreq = ft_getopt(cfg,'avgoverfreq','no');
cfg.clusterstatistic = 'maxsum';
cfg.clustertail = cfg.tail; 
% prepare neighbourhood structure
cfg_neighb.method = 'template'; % specify with which sensors other sensors can form clusters
cfg_neighb.layout = layout;
cfg.neighbours = ft_prepare_neighbours(cfg_neighb);
% minimum number of sensor neighbours in order for the cluster to be 
% accepted (see https://mailman.science.ru.nl/pipermail/fieldtrip/2011-February/003473.html)
cfg.minnbchan = 3; 

% Set up design matrix
nSubj = length(ftDataCell_cohZ_c1_vs_c2);
design = ones(2,2*nSubj);
design(1,nSubj+1:end) = 2;
design(2,:) = repmat(1:nSubj,1,2);
cfg.design = design;
cfg.ivar = 1; % independent variable
cfg.uvar = 2; % unit of observation variable

% Run the cluster-based permutation statistics
ftStat = ft_freqstatistics(cfg,ftDataCell_cohZ_c1_vs_c2{:},ftDataCell_coh_null{:});

end


function ftStat_aud = statTtest(ftDataCell_coh_c1,ftDataCell_coh_c2,cfg)
% Cluster permutation-based coherence statistics using the t-statistic

% Some additional tweaks to pretend the cohernece data are powerspectra, so
% ft_freqstatistics doesn't crash, see: 
% https://mailman.science.ru.nl/pipermail/fieldtrip/2018-October/038311.html
ftDataCell_coh_c1 = cellfun(@ftDataCohToFreq,ftDataCell_coh_c1,...
                                  'UniformOutput',false);
ftDataCell_coh_c2 = cellfun(@ftDataCohToFreq,ftDataCell_coh_c2,...
                                  'UniformOutput',false);

% Set up frequency statistic parameters
cfg.channel = ft_getopt(cfg, 'channel', 'all');
combinePlanar = any(~cellfun(@isempty,regexp(ftDataCell_coh_c1{1}.label,'.*\+.*','once')));
if strcmp(cfg.channel,'MEGMAG')
    layout = 'neuromag306mag.lay';
elseif strcmp(cfg.channel,'MEGGRAD') && ~combinePlanar
    layout = 'neuromag306planar.lay';
elseif strcmp(cfg.channel,'MEGGRAD') && combinePlanar
    layout = 'neuromag306cmb.lay';
end
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT'; 
cfg.numrandomization = ft_getopt(cfg,'numrandomization',5000);
cfg.correctm = 'cluster';
cfg.alpha = 0.05; 
cfg.tail = ft_getopt(cfg,'tail',1); 
if cfg.tail == 0,  cfg.correcttail = 'alpha'; end
cfg.clusteralpha = 0.05;
cfg.frequency = ft_getopt(cfg,'frequency',[2,8]); % Based on auditory and lip envelope coherence
cfg.avgoverfreq = ft_getopt(cfg,'avgoverfreq','no');
cfg.clusterstatistic = 'maxsum';
cfg.clustertail = cfg.tail; 
% prepare neighbourhood structure
cfg_neighb.method = 'template'; % specify with which sensors other sensors can form clusters
cfg_neighb.layout = layout;
cfg.neighbours = ft_prepare_neighbours(cfg_neighb);
% minimum number of sensor neighbours in order for the cluster to be 
% accepted (see https://mailman.science.ru.nl/pipermail/fieldtrip/2011-February/003473.html)
cfg.minnbchan = 3; 

% Set up design matrix
nSubj = length(ftDataCell_coh_c1);
design = ones(2,2*nSubj);
design(1,nSubj+1:end) = 2;
design(2,:) = repmat(1:nSubj,1,2);
cfg.design = design;
cfg.ivar = 1; % independent variable
cfg.uvar = 2; % unit of observation variable

% Run the cluster-based permutation statistics
ftStat_aud = ft_freqstatistics(cfg,ftDataCell_coh_c1{:},ftDataCell_coh_c2{:});

end


function ftData = ftDataCoherency2coherence(ftData)
ftData.cohspctrm = abs(ftData.cohspctrm);
end


function ftData = ftDataCohToFreq(ftData)
% Function to convert FT coherence data structures to freq data structures
ftData.label = ftData.labelcmb(:,1);
ftData.powspctrm = ftData.cohspctrm;
ftData = rmfield(ftData,{'labelcmb','cohspctrm'});
ftData.dimord = 'chan_freq';

end


function ftDataCohZ = ftDataCohZtransform(ftC1,ftC2)
% Wrapper function to do Z transformation of coherence on FT data structures
ftDataCohZ = ftC1;
ftDataCohZ.cohspctrm = cohZtransform(ftC1.cohspctrm,ftC2.cohspctrm,ftC1.dof,ftC2.dof);

end


function cohZ = cohZtransform(c1,c2,dfc1,dfc2)
% Function for computing the Z transformed coherence difference
% Snippet from ft_statfun_indepsamplesZcoh with some modifications
biasc1 = 1./(dfc1-2);
biasc2 = 1./(dfc2-2);
denomZ = sqrt(1./(dfc1-2) + 1./(dfc2-2));
cohZ = (atanh(c1)-biasc1-atanh(c2)+biasc2)./denomZ;

end