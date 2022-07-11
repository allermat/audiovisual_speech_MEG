function varargout = AVSM_setupdir(varargin)
% Set up folders for the experiment. 
% 
% USAGE: 
%   out = AVSM_setupdir();
%   out = AVSM_setupdir(dirID);
%   out = AVSM_setupdir(dirID,subID);
%
% INPUT:
%   dirID: directory ID string
%   subID: subject ID string
% OUTPUT: 
%   outPath: path of the required directory
%
% SIDE EFFECTS:
%   Some folders are created if not found. 
%   The required toolboxes are added to the path
%
% Copyright(C) Mate Aller 2020
% allermat@gmail.com

%% Parsing input. 
p = inputParser;

validDirIDs = {'analysis_megcoherence','analysis_megcoherence_group',...
    'analysis_megcoherence_sub','code_analysis','results','rawdata',...
    'rawdata_anat_sub','study_root','toolbox'};
validSubIDs = strcat('sub-',arrayfun(@(x) sprintf('%02d',x),1:14,...
    'UniformOutput',false));

addOptional(p,'dirID','',@(x) ismember(x,validDirIDs));
addOptional(p,'subID','',@(x) ismember(x,validSubIDs));

parse(p,varargin{:});

dirID = p.Results.dirID;
subID = p.Results.subID;

%% Setting up the basic directories if necessary. 
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[\w -]*','match');
setupID = setupID{1};

if strcmp(computer,'PCWIN64') || strcmp(computer,'PCWIN')
    [~,userID] = system('echo %username%');
    userID = regexp(userID,'[\w -]*','match');
    userID = userID{1};
elseif strcmp(computer,'GLNXA64') || strcmp(computer,'GLNX86')
    [~,userID] = system('id -n -u');
    userID = regexp(userID,'[\w -]*','match');
    userID = userID{1};
else
    error('Can''t find user name!');
end

% base folder depending on the setup
% You must add your setup to the list below and specify the path to the 
% repository
if strcmpi(setupID,'PC0220')
    baseDir = 'U:';
    mode = 'home';
elseif ~isempty(regexp(setupID,'^login','once'))
    baseDir = fullfile('/imaging/davis/users',userID);
    mode = 'analysis';
elseif ~isempty(regexp(setupID,'^node','once'))
    baseDir = fullfile('/imaging/davis/users',userID);
    mode = 'analysis';
else
    error('Unidentified setup!')
end

mypath.study_root = fullfile(baseDir,'Projects','AVSpeechMEG');
if ~exist(mypath.study_root,'dir')
    mkdir(mypath.study_root);
end

% mypath.analysis_behav = fullfile(mypath.study_root,'behavioural_analysis');
% if ~exist(mypath.analysis_behav,'dir') && strcmp(mode,'home')
%     mkdir(mypath.analysis_behav);
% end
% 
% mypath.analysis_eeg = fullfile(mypath.study_root,'EEG_analysis');
% if ~exist(mypath.analysis_eeg,'dir') && any(strcmp(mode,{'home','analysis'}))
%     mkdir(mypath.analysis_eeg);
% end

mypath.analysis_megcoherence = fullfile(mypath.study_root,'data','derivatives','megcoherence');
if ~exist(mypath.analysis_megcoherence,'dir') && any(strcmp(mode,{'home','analysis'}))
    mkdir(mypath.analysis_megcoherence);
end

mypath.code_analysis = fullfile(mypath.study_root,'code','analysis');
if any(strcmp(mode,{'home','analysis'}))
    if ~exist(mypath.code_analysis,'dir')
        mkdir(mypath.code_analysis);
    end
end

mypath.code_simulations = fullfile(mypath.study_root,'code','simulations');
if any(strcmp(mode,{'home','analysis'}))
    if ~exist(mypath.code_simulations,'dir')
        mkdir(mypath.code_simulations);
    end
end

mypath.results = fullfile(mypath.study_root,'results');
if any(strcmp(mode,{'home','analysis'}))
    if ~exist(mypath.results,'dir')
        mkdir(mypath.results);
    end
end

% mypath.data_behav = fullfile(mypath.study_root,expStage,'behavioural_data');
% if ~exist(mypath.data_behav,'dir') && any(strcmp(mode,{'home','presentation'}))
%     mkdir(mypath.data_behav);
% end

% mypath.data_meg = fullfile(mypath.study_root,'data','rawdata')
% mypath.data_meg = fullfile(mypath.study_root,expStage,'MEG_data');
% if ~exist(mypath.data_meg,'dir') && any(strcmp(mode,{'home','analysis'}))
%     mkdir(mypath.data_meg);
% end

mypath.rawdata = fullfile(mypath.study_root,'data','rawdata');
if ~exist(mypath.rawdata,'dir') && any(strcmp(mode,{'home','analysis'}))
    mkdir(mypath.rawdata);
end

mypath.toolbox = fullfile(mypath.study_root,'code','toolbox');
if ~exist(mypath.toolbox,'dir')
    mkdir(mypath.toolbox);
end

%% Adding folders to the path
% This is executed only when the function is called without input
if isempty(varargin)
    % Study root folder
    addpath(fullfile(mypath.study_root));
    % All folders in toolbox
    d = dir(mypath.toolbox);
    dirNames = {d.name}';
    dirNames = dirNames([d.isdir]);
    dirNames = dirNames(~ismember(dirNames,{'.','..'}));
    dirPaths = fullfile(mypath.toolbox,dirNames);
    addpath(mypath.toolbox);
    addpath(dirPaths{:});
    
    % Initialize fieldtrip if it is added to the path
    if ismember('fieldtrip',dirNames)
        ft_defaults;
    end
    
    if any(strcmp(mode,{'home','analysis'}))
        % analysis scripts
        addpath(mypath.code_analysis);
        % simulations
        addpath(mypath.code_simulations);
        % results
        addpath(mypath.results);
    end
end

%% Returning the required path
% If not found either an error is thrown or the folder is created. 

if strcmp(dirID,'analysis_megcoherence')
    outPath = mypath.analysis_megcoherence;
elseif strcmp(dirID,'analysis_megcoherence_group')
    outPath = fullfile(mypath.analysis_megcoherence,'group');
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'analysis_megcoherence_sub')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.analysis_megcoherence,subID);
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'code_analysis')
    outPath = mypath.code_analysis;
% elseif strcmp(dirID,'data_behav')
%     outPath = mypath.data_behav;
% elseif strcmp(dirID,'data_meg')
%     outPath = mypath.data_meg;
elseif strcmp(dirID,'rawdata_anat_sub')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.rawdata,subID,'anat');
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
% elseif strcmp(dirID,'data_behav_sub')
%     if strcmp(subID,'')
%         error('Subject ID must be specified!');
%     end
%     outPath = fullfile(mypath.data_behav,subID);
%     if ~exist(outPath,'dir')
%         mkdir(outPath);
%     end
% elseif strcmp(dirID,'data_meg_sub')
%     if strcmp(subID,'')
%         error('Subject ID must be specified!');
%     end
%     outPath = fullfile(mypath.data_meg,subID);
%     if ~exist(outPath,'dir')
%         mkdir(outPath);
%     end
elseif strcmp(dirID,'results')
    outPath = mypath.results;
elseif strcmp(dirID,'study_root')
    outPath = mypath.study_root;
elseif strcmp(dirID,'toolbox')
    outPath = mypath.toolbox;
else
    outPath = '';
end

if ~isempty(outPath)
    varargout{1} = outPath;
else
    fprintf('\nSetup ID: %s\nUser ID: %s\n\n',setupID,userID);
end

end