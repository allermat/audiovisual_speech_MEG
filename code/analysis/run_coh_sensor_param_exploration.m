clearvars;
% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};
fprintf(setupID);
if ~isempty(regexp(setupID,'^node.+','once'))
    runMode = 'node';
elseif ~isempty(regexp(setupID,'^login.+','once'))
    runMode = 'login';
    dbstop if error;
else
    runMode = 'other';
    dbstop if error;
end

% Opening parallel pool. 
% if there is no parallel pool running, open one. 
currPool = gcp('nocreate');
if isempty(currPool)
    switch runMode
        case 'node'
            parpool('local',16);
        case 'login'
            parpool('local',4);
        case 'other'
            parpool('local');
    end
end

% Compute coherence for each subject
subjList = strcat('sub-',arrayfun(@(x) sprintf('%02d',x),1:14,'UniformOutput',false));
analysisDetails = logical(fullfact([2,2])-1);
ftStatCell = cell(size(analysisDetails,1),2);
for i = 1:size(analysisDetails,1)
    parfor iSubj = 1:numel(subjList)
        coh_sensor_aud_vis(subjList{iSubj},'cohAnalyses','basicAudVis',...
            'doDerangement',analysisDetails(i,1),...
            'doSquaredCoh',analysisDetails(i,2));
    end
    % Compute group level coherence statistics on basic audio and visual
    % coherence effects
    ftStatCell(i,:) = coh_sensor_group_stats('basicAudVis');
end
results = table(analysisDetails(:,1),analysisDetails(:,2),ftStatCell(:,1),...
                ftStatCell(:,2),'VariableNames',{'deranged','squared',...
                'basicAud','basicVis'});
destDir = AVSM_setupdir('analysis_megcoherence_group');
save(fullfile(destDir,'coh_sensor_stats_param_exploration.mat'),'results','-v7.3');

%% Optionally, plot figures
plotFigures = false;
if plotFigures
    destDir = AVSM_setupdir('analysis_megcoherence_group');
    load(fullfile(destDir,'coh_sensor_stats_param_exploration.mat'),'results');
    for i = 1:size(results,1)
        cfg = [];
        cfg.alpha     = 0.05;
        cfg.parameter = 'stat';
        cfg.layout    = 'neuromag306mag.lay';
        cfg.subplotsize = [1,5];
        cfg.zlim = 'maxabs';
        ft_clusterplot(cfg, results.basicAud{i});
        colorbar;
        set(gcf, 'Position', [0 0 1000 1500], 'PaperPositionMode', 'auto');
        
        cfg = [];
        cfg.alpha     = 0.05;
        cfg.parameter = 'stat';
        cfg.layout    = 'neuromag306mag.lay';
        cfg.subplotsize = [1,5];
        cfg.zlim = 'maxabs';
        ft_clusterplot(cfg, results.basicVis{i});
        colorbar;
        set(gcf, 'Position', [0 0 1000 1500], 'PaperPositionMode', 'auto');
    end
end
    