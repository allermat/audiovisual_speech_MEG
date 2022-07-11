%% Import data from spreadsheet
df_path = fullfile(AVSM_setupdir('study_root'),'data','derivatives',...
                   'behavioural','group');
% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 6);

% Specify sheet and range
opts.Sheet = "Average results";
opts.DataRange = "A2:F15";

% Specify column names and types
opts.VariableNames = ["VarName1", "AO20", "AO70", "AV20", "AV70", "VO"];
opts.SelectedVariableNames = ["VarName1", "AO20", "AO70", "AV20", "AV70", "VO"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");

% Import the data
wordReportOrig = readtable(fullfile(df_path,'AVSpeechMEG_WordReport_allsubj_osf.xls'),...
                           opts,"UseExcel",false);

%% Do some tidying
% The first column contains the subject IDs as the data were collected. 
% Some participants were rejected before data analysis so I recoded the
% subject IDs such that they go from 1-14. 
temp = wordReportOrig.Properties.VariableNames;
temp{1} = 'subID';
wordReportOrig.Properties.VariableNames = temp;
temp = strcat('sub-',arrayfun(@(x) sprintf('%02d',x),(1:14)','UniformOutput',false));
wordReportOrig.subID = temp;

%% Transform the table to a more widly compatible format
% This way word report accuracy becomes a single variable (column) and
% other categorical variables are added (stim_modality, acc_clarity)
tempSubj = repmat(wordReportOrig.subID,5,1);
tempAcc = wordReportOrig{:,2:6};
tempAcc = tempAcc(:);
tempModality = repmat({'AO','AO','AV','AV','VO'},14,1);
tempModality = categorical(tempModality(:));
tempClarity = repmat({'low','high','low','high','none'},14,1);
tempClarity = categorical(tempClarity(:));

wordReport = table(tempSubj,tempModality,tempClarity,tempAcc,'VariableNames',...
                   {'subID','stim_modality','acc_clarity','accuracy'});
               
% Save data in various formats for later use
save(fullfile(df_path,'AVSpeechMEGWordReportAllSubj.mat'),...
     'wordReportOrig','wordReport')
writetable(wordReport,fullfile(df_path,'word_report.csv'));
word_report_data = wordReportOrig{:,2:end};
save(fullfile(df_path,'AVSpeechMEGWordReportAllSubj_python.mat'),...
     'word_report_data')