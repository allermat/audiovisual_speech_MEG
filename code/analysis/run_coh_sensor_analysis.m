%% Prepare parallel pool
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
compIndividual = true;
doGroupStats = true;

%% Compute coherence for each subject
if compIndividual
    subjList = strcat('sub-',arrayfun(@(x) sprintf('%02d',x),1:14,'UniformOutput',false));
    parfor i = 1:numel(subjList)
        coh_sensor_aud_vis(subjList{i},'doDerangement',true);
    end
end
%% Compute group statistics
if doGroupStats
    % Compute group level coherence statistics on basic audio and visual
    % coherence effects
    ftStatCell = coh_sensor_group_stats('basicAudVis');
    
    % Plot basic audio and visual coherence effects
    cmap = crameri('vik');
    h = ceil(size(cmap,1)/2);
    cfg = [];
    cfg.alpha     = 0.05;
    cfg.parameter = 'stat';
    cfg.layout    = 'neuromag306mag.lay';
    cfg.subplotsize = [1,5];
    cfg.zlim = [tinv(0.95,13),tinv(0.9998,13)];
    cfg.style = 'straight';
    ft_clusterplot(cfg, ftStatCell{1});
    colormap(flipud(cmap(1:h,:)));
    colorbar;
    set(gcf, 'Position', [0 0 1000 1500], 'PaperPositionMode', 'auto');
    
    cfg = [];
    cfg.alpha     = 0.05;
    cfg.parameter = 'stat';
    cfg.layout    = 'neuromag306mag.lay';
    cfg.subplotsize = [1,5];
    cfg.zlim = [tinv(0.95,13),tinv(0.9998,13)];
    cfg.style = 'straight';
    ft_clusterplot(cfg, ftStatCell{2});
    colormap(cmap(h+1:end,:));
    colorbar;
    set(gcf, 'Position', [0 0 1000 1500], 'PaperPositionMode', 'auto');
    
end
