function ranovaTable = behav_word_report_rmanova(wordReportOrig)
% Perform repeated measures ANOVA on word report accuracy data

% Copyright(C) 2016, Mate Aller
%  allermat@gmail.com

% Input data structure: 
% Columns: responses ordered according to the conditions of the within
% subjects design. Basically these are the repeated measurements.  
% Rows: subjects. 
wordReportOrig = wordReportOrig(:,2:5);
wordReportOrig.Properties.VariableNames = {'AO_20','AO_70','AV_20','AV_70'};

% First we have to fit a repeated measures model and then we perform the
% repeated measures anova on the fitted model. To fit the model we have to
% specify the within design. This is basically a table with the full 
% factorial expansion of the two within factors in our case. The number of 
% rows has to match the number of response variables in the data. 
withinDesign = table(...
    categorical({'AO';'AO';'AV';'AV'}),...
    categorical({'20';'70';'20';'70'}),...
    'VariableNames',{'modality','clarity'});
% To fit the model we feed the available data, the between subjects model
% and the within subject design into the function. 
% The input data:
% should contain the response (within subject) and explanatory (between 
% subjects) variables. In our case we don't have between subjects variables
% (all subjects are equal), hence we just include the responses here. 
% The between subjects model:
% is given using the Wilkinson notation (you can find more on this in the 
% documentation). Basically it is the classical Y ~ X1 + X2 notation 
% (where Y: response, X1,X2 explanatory variables). In this particular case
% because of the repeated measures design, we have six response variables 
% (corresponding to the full factorial expansion of the within subjects 
% factors). Also because we don't have any between subjects factors we just
% use the intercept term as explanatory (dummy) variable. 
% Hence our between subjects model looks like this: Y1,Y2,Y3,Y4,Y5,Y6 ~ 1 
% (where Y1-Y6 response variables). What you see in the input of the 
% function is equivalent to this just using the actual names of the 
% variables and a simplified way of specifying a list of variables 
% (again, see the documentation for details). 
% The within subjects design: 
% we're just including the prepaired table. 
rm = fitrm(wordReportOrig,'AO_20-AV_70 ~ 1','WithinDesign',withinDesign);
% The output of the function is a RepeatedMeasuresModel object. It computes
% a bounch of things, have a look at the documentation for the details. 

% Now we're at the point of actually performing the repeated measures
% ANOVA. Just call the ranova function with the repeated measures model 
% object and specify the within subjects model (again using the Wilkinson
% notation). In our case we are interested in the main effects and the 
% interaction of the within subjects factors, so I defined this there. 
ranovaTable = ranova(rm,'WithinModel','modality*clarity');

end
