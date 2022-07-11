function [s, cfg] = ft_statfun_depsamplesDiff(cfg, dat, design)

% FT_STATFUN_DIFF computes the difference of the mean in two conditions. Although it
% can be used for statistical testing, it is not very usefull since it will have
% rather limited sensitivity.
% 
% The purpose of this function is to show you an example on how to write a statfun
% that expresses the difference in the data between two conditions. You can use such
% a function with the statistical framework in FieldTrip to perform a simple (or more
% complex) permutation test, without having to worry about the representation of the
% data.
%
% See also FT_STATFUN_MEAN for another example function

% Copyright (C) 2006, Robert Oostenveld 
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% perform some checks on the configuration
if ~isfield(cfg,'uvar') || isempty(cfg.uvar)
  ft_error('uvar must be specified for dependent samples statistics');
end

% perform some checks on the design
sel1 = find(design(cfg.ivar,:)==1);
sel2 = find(design(cfg.ivar,:)==2);
n1  = length(sel1);
n2  = length(sel2);
if (n1+n2)<size(design,2) || (n1~=n2)
  ft_error('Invalid specification of the design array.');
end
nunits = length(design(cfg.uvar, sel1));
df = nunits - 1;
if nunits<2
  ft_error('The data must contain at least two units (usually subjects).')
end
if (nunits*2)~=(n1+n2)
  ft_error('Invalid specification of the design array.');
end

% First compute the differences within units of observation
% store the positions of the 1-labels and the 2-labels in a nunits-by-2 array
poslabelsperunit = zeros(nunits,2);
poslabel1        = find(design(cfg.ivar,:)==1);
poslabel2        = find(design(cfg.ivar,:)==2);
[dum,i]          = sort(design(cfg.uvar,poslabel1), 'ascend');
poslabelsperunit(:,1) = poslabel1(i);
[dum,i]          = sort(design(cfg.uvar,poslabel2), 'ascend');
poslabelsperunit(:,2) = poslabel2(i);

% calculate the differences between the conditions
diffmat = dat(:,poslabelsperunit(:,1)) - dat(:,poslabelsperunit(:,2));

% Then compute the mean over the differences
s.stat = nanmean(diffmat,2);

end

% the stat field is used in STATISTICS_MONTECARLO to make the
% randomization distribution, but you can also return other fields
% which will be passed on to the command line in the end.

