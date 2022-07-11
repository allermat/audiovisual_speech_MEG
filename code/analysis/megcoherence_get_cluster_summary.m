function cluster_summary = megcoherence_get_cluster_summary(fName)

sourceDir = AVSM_setupdir('analysis_megcoherence_group');
dataFile = fullfile(sourceDir,fName);
load(dataFile);
% Keep in mind that all indices are from python where indexing starts from
% 0, so here +1 must be added to be correct
nGoodCluster = numel(good_cluster_inds);
% The first half of vertices belong to the left hemisphere, the second to
% the right
idxLeftHemi = 1:size(T_obs,2)/2;
idxRightHemi = (size(T_obs,2)/2)+1:size(T_obs,2);
if iscolumn(good_cluster_p), good_cluster_p = good_cluster_p'; end
for i = nGoodCluster:-1:1
    temp = clusters{good_cluster_inds(i)+1,2};
    if all(ismember(temp,idxLeftHemi))
        hemi{i} = 'lh';
    elseif all(ismember(temp,idxRightHemi))
        hemi{i} = 'rh';
    else
        error('Cluster belongs to more than one hemisphere')
    end
    nvtx(i) = numel(temp);
    t_sum(i) = sum(T_obs(temp+1));
end

cluster_summary = table(hemi',nvtx',t_sum',good_cluster_p','VariableNames',...
                        {'hemisphere','n_vertices','t_sum','p_val'});
                    