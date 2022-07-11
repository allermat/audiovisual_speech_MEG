%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Trans head positions 
%
%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Olaf Hauk (March 2009), modified by Lucy MacGregor and Heidi
% Solberg Okland
%
% http://imaging.mrc-cbu.cam.ac.uk/meg/InterpolateData
%
% Script interpolates data from multiple files to a specified sensor array.
% It is used to normalise data across subjects who will have different head
% positions, thus allowing for data analysis at sensor level.
%
% 1. Specify files to be interpolated in cell array "fiff_files{}" 
% 2. Specify sensor array to be used in "trans_file"
% 3. Apply maxfilter with -trans option to files specified in 1. to
% interpolate data to sensor array specified in 2.
% 
% Original file names will be appended by whatever is specified in "t_ext"
% For each trans-ed file, a log-file with maxfilter output will be created ("*.log")
%
% NB. "trans_file" can be the .fif file determined as having sensor array
% closest to the average of your set of files (see "average_sensorarray.m")

% Set up datapath and files
dataPath = '/imaging/hs01/AV_speech/M+EEG'; cd(dataPath)

% Specify subjects
subjPaths = {};

% Generic filename for the input .fif-file
file_in = ['concatenated_icaed_raw.fif'];
t_ext = '_trans';   % suffix for filenames of trans-ed files

% Specify the fiff file containing the desired headposition - in this case 
% we want the subject with the headposition that is closest to the average 
% (calculated in headpos_average.m)
trans_file = [dataPath '' file_in]; % .fif file with sensor geometry to which other files will be interpolated

% Get all the files..
for ss = 1:length(subjPaths),  
    fiff_files{ss} = strjoin(strcat(dataPath,subjPaths(ss),file_in));
end;

% Check if trans file exists..
if ~exist(trans_file, 'file'),
    fprintf(1, 'Cannot proceed without trans_file %s\n', trans_file);
    return;
end;

trans_data = fiff_setup_read_raw(trans_file); % read in raw data
nr_files = length(fiff_files);

% For each input file, read in, apply maxfilter -trans, write the output fiff file and log file
for ff = 1:nr_files, 
    
    fprintf('Applying maxfilter -trans on subject %d of %d...\n', ff, nr_files)
    
    in_file = fiff_files{ff};
    [thispath,thisfile,thisext] = fileparts(in_file);
    out_file = fullfile(thispath, [thisfile t_ext thisext]);
    log_file = fullfile(thispath, [thisfile t_ext '.log']);
    
    % Go!
    eval(sprintf('!/neuro/bin/util/maxfilter -f %s -o %s -trans %s -v -frame head -origin 0 0 45 -force -ctc /neuro/databases/ctc/ct_sparse.fif -cal /neuro/databases/sss/sss_cal.dat | tee %s', in_file, out_file, trans_file, log_file));

end;
