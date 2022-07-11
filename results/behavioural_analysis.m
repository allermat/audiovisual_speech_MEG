%% Load data
load(fullfile(AVSM_setupdir('study_root'),'data','derivatives','behavioural',...
              'group','AVSpeechMEGWordReportAllSubj.mat'),'wordReport',...
              'wordReportOrig');

%% Explore plotting options
clear g

g(1,1)=gramm('x',wordReport.stim_modality,'y',wordReport.accuracy,...
             'color',wordReport.acc_clarity);
g(1,1).set_order_options('x',{'VO','AO','AV'},'color',{'none','high','low'});
g(1,1).set_color_options('map','brewer_dark');
g(1,2)=copy(g(1));
g(1,3)=copy(g(1));
g(2,1)=copy(g(1));
g(2,2)=copy(g(1));
g(2,3)=copy(g(1));

%Raw data as scatter plot
g(1,1).geom_point();
g(1,1).set_title('geom_point()');

%Jittered scatter plot
g(1,2).geom_jitter('width',0.4,'height',0);
g(1,2).set_title('geom_jitter()');

%Averages with confidence interval
g(1,3).stat_summary('type','ci','geom',{'line','point','errorbar'});
g(1,3).set_title('stat_summary()');
%Boxplots
g(2,1).stat_boxplot();
g(2,1).set_title('stat_boxplot()');

%Violin plots
g(2,2).stat_violin('fill','transparent');
g(2,2).set_title('stat_violin()');

%Jittered scatter plot
g(2,3).geom_jitter('width',0.4,'height',0);
%Averages with confidence interval
g(2,3).stat_summary('type','ci','geom',{'line','black_errorbar'});
g(2,3).set_title('geom_jitter() + stat_summary()');

%These functions can be called on arrays of gramm objects
g.set_names('x','Modality','y','Accuracy','color','Clarity');
g.set_title('Word report accuracies');

figure('Position',[100 100 800 550]);
g.draw();

%% AO70 vs AV20
clear g1

idx1 = wordReport.stim_modality == 'AO' & wordReport.acc_clarity == 'high';
idx2 = wordReport.stim_modality == 'AV' & wordReport.acc_clarity == 'low';
col = wordReport.accuracy(idx2)-wordReport.accuracy(idx1);
g1(1,1)=gramm('x',wordReport.stim_modality(idx1 | idx2),'y',wordReport.accuracy(idx1 | idx2),...
              'color',repmat(col,2,1));
g1(1,1).no_legend();
%Raw data as scatter plot
g1(1,1).geom_point();
g1(1,1).geom_line();
g1(1,1).set_title('geom_point()');

g1.set_names('x','Modality','y','Accuracy');
g1.set_title('Word report accuracies');

figure();
g1.draw();

% Simple matlab plot
x = repmat([1,2],14,1);
x = x(:);
y = cat(1,wordReport.accuracy(idx1),wordReport.accuracy(idx2));
hFig = figure();
colormap(hFig,'redblue');
scatter(x,y,[],repmat(col,2,1),'filled','LineWidth',2,'Marker','o');
xlim([0,3]);
set(gca,'xTick',[1,2],'xTickLabels',unique(wordReport.stim_modality(idx1 | idx2)));

%% Working figure
clear g2

idx1 = wordReport.stim_modality == 'AO' & wordReport.acc_clarity == 'high';
idx2 = wordReport.stim_modality == 'AV' & wordReport.acc_clarity == 'low';
col = wordReport.accuracy(idx2)-wordReport.accuracy(idx1);

% All conditions
g2(1,1)=gramm('x',wordReport.stim_modality,'y',wordReport.accuracy,...
             'color',wordReport.acc_clarity);
g2(1,1).set_order_options('x',{'VO','AO','AV'},'color',{'none','high','low'});
g2(1,1).set_color_options('map','brewer_dark');
% Jittered scatter plot
g2(1,1).geom_jitter('width',0.4,'height',0);
% Averages with confidence interval
g2(1,1).stat_summary('type','sem','geom',{'line','black_errorbar'});
g2(1,1).set_title('All conditions');
g2(1,1).set_layout_options('legend_position',[0.1,0.5,0.3,0.4]);
g2(1,2)=gramm('x',wordReport.stim_modality(idx1 | idx2),'y',wordReport.accuracy(idx1 | idx2),...
              'color',repmat(col,2,1));

% AO70 vs AV20
g2(1,2).no_legend();
g2(1,2).geom_point();
g2(1,2).geom_line();
g2(1,2).set_title('AO_High vs AV_Low');
g2.set_names('x','Modality','y','Accuracy','color','Clarity');
g2.set_title('Word report accuracies');

figure('Position',[100 100 800 400]);
g2.draw();

%% Some descriptive statistics for AV, AO, and VOclo
temp_data = wordReport(ismember(wordReport.stim_modality,{'AO','AV'}),:);
descr = join(varfun(@mean,temp_data,'InputVariables','accuracy',...
                    'GroupingVariables',{'stim_modality','acc_clarity'}),...
             varfun(@std,temp_data,'InputVariables','accuracy',...
                    'GroupingVariables',{'stim_modality','acc_clarity'}));
temp_data = wordReport(ismember(wordReport.stim_modality,{'VO'}),:);
temp = join(varfun(@mean,temp_data,'InputVariables','accuracy',...
                    'GroupingVariables',{'stim_modality','acc_clarity'}),...
             varfun(@std,temp_data,'InputVariables','accuracy',...
                    'GroupingVariables',{'stim_modality','acc_clarity'}));
descr = cat(1,descr,temp);
descr.sem_accuracy = descr.std_accuracy./sqrt(descr.GroupCount);
%% Repeated measures ANOVA on word report data
ranovaTable = behav_word_report_rmanova(wordReportOrig) %#ok<NOPTS>

vis_speech_eff = join(varfun(@mean,wordReport(ismember(wordReport.stim_modality,{'AO','AV'}),:),...
                             'InputVariables',{'accuracy'},'GroupingVariables',{'stim_modality'}),...
                      varfun(@std,wordReport(ismember(wordReport.stim_modality,{'AO','AV'}),:),...
                             'InputVariables',{'accuracy'},'GroupingVariables',{'stim_modality'}));
vis_speech_eff.sem_accuracy = vis_speech_eff.std_accuracy./sqrt(vis_speech_eff.GroupCount/2);

acc_clarity_eff = join(varfun(@mean,wordReport(ismember(wordReport.stim_modality,{'AO','AV'}),:),...
                             'InputVariables',{'accuracy'},'GroupingVariables',{'acc_clarity'}),...
                       varfun(@std,wordReport(ismember(wordReport.stim_modality,{'AO','AV'}),:),...
                             'InputVariables',{'accuracy'},'GroupingVariables',{'acc_clarity'}));
acc_clarity_eff.sem_accuracy = acc_clarity_eff.std_accuracy./sqrt(acc_clarity_eff.GroupCount/2);

% Compute simple effects as the interaction is significant
res = cell(4,4);
% AV vs AO in high clarity
[res{1,:}] = ttest(wordReportOrig.AV70,wordReportOrig.AO70);
% AV vs AO in low clarity
[res{2,:}] = ttest(wordReportOrig.AV20,wordReportOrig.AO20);
% high vs. low clarity in AO
[res{3,:}] = ttest(wordReportOrig.AO70,wordReportOrig.AO20);
% high vs. low clarity in AV
[res{4,:}] = ttest(wordReportOrig.AV70,wordReportOrig.AV20);
% Organise results in a neat table
temp = cellfun(@(x) struct2cell(x)',res(:,4),'UniformOutput',false);
temp = cat(1, temp{:});
simple_effects = cell2table(cat(2,res(:,1:3),temp),'VariableNames',...
                            {'h','p','ci','tstat','df','sd'},'RowNames', ...
                            {'AV_vs_AO_in_high','AV_vs_AO_in_low',...
                             'high_vs_low_in_AO','high_vs_low_in_AV'}) %#ok<NOPTS>
                         
%% Bayes-factors for AVlow vs AOhigh
X = wordReport.accuracy(ismember(wordReport.stim_modality,'AV') & ...
                        ismember(wordReport.acc_clarity,'low'));
Y = wordReport.accuracy(ismember(wordReport.stim_modality,'AO') & ...
                        ismember(wordReport.acc_clarity,'high'));
[bf10,p,CI,stats] = bf.ttest(X,Y);
bf01 = 1/bf10;
bayesian_ttest_AVlow_AOhigh = table(bf01,bf10,p,CI',stats.tstat,stats.df,stats.sd,...
                                 'VariableNames',{'bf01','bf10','p','ci','tstat','df','sd'});

%% Correlation between (AV20-AO70) and VO
X = cat(2,wordReportOrig.VO,wordReportOrig.AV20-wordReportOrig.AO70);
figure('Name','Correlation AV_low-AO_high vs. VO');
[R,p_val] = corrplot(X,'varNames',{'VO','AV-AO'},'testR','on'); 
temp = struct();
temp.var1 = {'AV_low-AO_high'};
temp.var2 = {'VO'};
temp.R = R(2,1);
temp.p_val = p_val(2,1);
vis_pref_index_corr = struct2table(temp);
% Correlation between (AV-AO) and VO as per request of reviewer #1
X = cat(2,wordReportOrig.VO,...
          mean([wordReportOrig.AV20,wordReportOrig.AV70],2) - ...
          mean([wordReportOrig.AO20,wordReportOrig.AO70],2));
figure('Name','Correlation AV-AO vs. VO');
[R,p_val] = corrplot(X,'varNames',{'VO','AV-AO'},'testR','on'); 
temp = struct();
temp.var1 = {'AV-AO'};
temp.var2 = {'VO'};
temp.R = R(2,1);
temp.p_val = p_val(2,1);
vis_pref_index_corr = cat(1,vis_pref_index_corr,struct2table(temp));

%% Correlation between (AV20-AO70) and VO with marginal densities
clear g
figure('Position',[100 100 550 550]);

ks_bandwidth = 0.07;
% Tha basic script comes from the gramm examples. 
%Create x data histogram on top
g(1,1)=gramm('x',wordReportOrig.VO);
g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.1 0.02],...
    'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
g(1,1).set_names('x','');
g(1,1).stat_density('bandwidth',ks_bandwidth); %ksdensity
% g(1,1).stat_bin('geom','line','fill','all','nbins',15); %histogram
g(1,1).axe_property('XTickLabel',''); % We deactivate tht ticks

%Create a scatter plot
g(2,1)=gramm('x',wordReportOrig.VO,...
             'y',wordReportOrig.AV20-wordReportOrig.AO70);
g(2,1).set_names('x','VO','y','AVlow-AOhigh');
g(2,1).geom_point(); % Scatter plot
g(2,1).stat_glm('disp_fit',false); % Fit line
g(2,1).set_point_options('base_size',6);
g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
    'legend_pos',[0.83 0.75 0.2 0.2],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);

%Create y data histogram on the right
g(3,1)=gramm('x',wordReportOrig.AV20-wordReportOrig.AO70);
g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
    'legend',false,...
    'margin_height',[0.1 0.02],...
    'margin_width',[0.02 0.05],...
    'redraw',false);
g(3,1).set_names('x','');
g(3,1).stat_density('bandwidth',ks_bandwidth); %ksdensity
% g(3,1).stat_bin('geom','line','fill','all','nbins',15); %histogram
g(3,1).coord_flip();
g(3,1).axe_property('XTickLabel','');

%Set global axe properties
g.axe_property('TickDir','out');
g.set_title('Correlation across subjects');
g.set_color_options('map','d3_10');
g.draw();

