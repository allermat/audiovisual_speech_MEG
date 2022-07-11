% Save envelopes and meg data separately for MNE-python analysis
clearvars;
subjList = dir(fullfile(AVSM_setupdir('analysis_megcoherence'),'sub-*'));
subjList = {subjList.name}';

for iSub = 1:numel(subjList)
    saveDir = cd(AVSM_setupdir('analysis_megcoherence_sub',subjList{iSub}));
    load('trials_with_env_MEG.mat')
    
    envelopes = cellfun(@(x) shiftdim(x(end-1:end,:),-1),...
                        trials_with_env_MEG.trial,'UniformOutput',false);
    envelopes = cat(1,envelopes{:});
    save('envelopes.mat','envelopes');
    
    cfg = [];
    cfg.channel = 'MEG';
    trials_no_env_MEG = ft_selectdata(cfg,trials_with_env_MEG);
    save('trials_no_env_MEG.mat','trials_no_env_MEG')
    
    cd(saveDir);
    
end


