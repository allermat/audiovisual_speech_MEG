function executeFunGroup(fun,matchStr,varargin)
% Runs mvpa functions on the group of subjects

% Parsing input, checking matlab
p = inputParser;


addRequired(p,'fun',@(x) validateattributes(x,{'function_handle'},{'scalar'}));
addRequired(p,'fileMatchStr',@ischar);
addParameter(p,'funInput',{},@(x) validateattributes(x,{'cell'},{'vector','row'}));
addParameter(p,'nFilesExpected',[],@(x) validateattributes(x,{'numeric'},{'scalar'}));

parse(p,fun,matchStr,varargin{:});

fun = p.Results.fun;
fileMatchStr = p.Results.fileMatchStr;
funInput = p.Results.funInput;
nFilesExpected = p.Results.nFilesExpected;

subjList = dir(fullfile(AVSM_setupdir('analysis_megcoherence'),'sub-*'));
subjList = {subjList.name}';

for iSub = 1:numel(subjList)
    fun(subjList{iSub},fileMatchStr,funInput{:})
end

end

function filePathList = collectFiles(subjList,expStage,dirID,fileMatchStr,trMethod,nFilesExpected,subFolder)

filePathList = {};
index = 1;
for i = 1:size(subjList,1)
    saveDf = cd(fullfile(DEC_2_setupdir(expStage,dirID,subjList{i}),trMethod,subFolder));
    fileList = cellstr(ls);
    matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));
    if sum(matchID) == 0
        warning('No file, skipping subject %s! ',subjList{i});
        cd(saveDf);
        continue;
    elseif sum(matchID) > nFilesExpected
        warning('More files than needed, skipping subject %s! ',subjList{i});
        cd(saveDf);
        continue;
    else
        fileName = fileList(matchID);
        if iscolumn(fileName)
            fileName = fileName';
        end
    end
    temp = cellfun(@fullfile,repmat({pwd},size(fileName)),fileName,'UniformOutput',false);
    filePathList = cat(1,filePathList,temp);
    index = index + 1;
    cd(saveDf);
end

if isrow(filePathList)
    filePathList = filePathList';
end

end