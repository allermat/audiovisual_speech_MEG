%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute average sensor array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Olaf Hauk (March 2009), modified by Lucy MacGregor and Heidi
% Solberg Okland
%
% http://imaging.mrc-cbu.cam.ac.uk/meg/StandardSensorArray
%
% Script to calculate average sensor transformation from a list of input
% files and identify the file which is closest to the average (according to
% both translation and rotation parameters).
%
% NB1. Recommended to use the file closest to the average for subsequent
% transformation of other files (see "trans_171114.m")

% addpath( pwd ) % for local functions

dataPath = '';
cd(dataPath)

subjPaths = {};

% Generic filename for the input .fif-file
file_in = ['concatenated_icaed_raw.fif'];

% Loop over subjects to extract headposition info and make individual averages
nr_files = length(subjPaths);
trans_gm = 0;
for f = 1:nr_files, % read all files and average all that's relevant
    fprintf(1, 'Processing %s...\n', subjPaths{f});
    
    % Specify current .fif-file
    fiff = strjoin(strcat(dataPath,subjPaths(f),file_in));
    
    % Read in data
    raw = fiff_setup_read_raw(fiff); % read in raw data
    %[data,times] = fiff_read_raw_segment(raw); % "data" is a matrix of channels x amplitude readings
    
    % Pick out information about headposition
    trans_gm = trans_gm + raw.info.dev_head_t.trans;  % average sensor transformation matrix
    all_trans(:,:,f) = raw.info.dev_head_t.trans;  % remember individual transformations

end;

% Average sensor transformation matrices
trans_gm = trans_gm/nr_files;

% Insert this average info into a copied data variable
data_new = raw;
data_new.info.dev_head_t.trans = trans_gm;

% Find input file that is closest to average...
for f = 1:nr_files,
    diff_trans_t = all_trans(1:3,4,f) - trans_gm(1:3,4);        % translation difference
    diff_trans_r = all_trans(1:3,1:3,f) - trans_gm(1:3,1:3);    % rotation difference
    norm_t(f) = norm( diff_trans_t );
    norm_r(f) = norm( diff_trans_r );
end;

[y_t, i_t] = sort( norm_t );    % sort according to translation difference
[y_r, i_r] = sort( norm_r );    % sort according to rotation difference

for f = 1:nr_files, % combine the two sorted lists
    rank_t(f) = find( (i_t-f)==0 ); % rank for translation
    rank_r(f) = find( (i_r-f)==0 ); % rank for rotation
    avg_rank(f) = (rank_t(f) + rank_r(f))/2;
    [min_val, min_file] = min( avg_rank );  % best average rank
end;

winner = subjPaths{min_file}; % filename of file with best average ran

fprintf('\nAnd the winner is... %s !\nThis subject''s headposition is the one that comes closest to the average.\n', winner);
