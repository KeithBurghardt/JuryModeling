function [y, t] = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h, boundaries)

% BFINITEEXPFIT fits a EVL model to the data from a sampled lambda distribution
% Source: http://tuvalu.santafe.edu/~aaronc/powerlaws/bins/
%
% Version 1.0 (2015)
% Copyright (C) 2015 Keith Burghardt (University of Maryland, College Park)
% Distributed under GNU GPL v3.0
% http://www.gnu.org/copyleft/gpl.html
% BEXPFIT comes with ABSOLUTELY NO WARRANTY
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------

NumTrials = 100000;

t = reshape(boundaries,numel(boundaries),1);
N = numel(h(:,1)) - 1;
y = zeros(numel(LambdaRangeGuilty),numel(LambdaRangeInnocent),N+1,numel(t));

FractVotes = sum(h,2)./sum(reshape(h,numel(h),1));
%NumGuiltyVotes = (0:N).*NumVotes;
%NumInnocentVotes = (N:-1:0).*NumVotes;
count_i = 0;
for i_guilty = LambdaRangeGuilty
    count_j = 0;
    count_i = count_i + 1; 
    for i_innocent = LambdaRangeInnocent
        count_j = count_j + 1; 
        for j = 0:N
            x = max([exprnd(i_guilty,j,NumTrials);exprnd(i_innocent,N-j,NumTrials)]); 
            h2 = histc(x, t);
            h2(end) = [];
            % Compute distance using KS statistic
            temp = cumsum(h2(end:-1:1));
            cx = [1 - temp(end:-1:1)./NumTrials 1];
            y(count_i,count_j,j+1,:) = FractVotes(j+1).*cx;%ksdensity(max([exprnd(i_guilty,j,NumTrials);exprnd(i_innocent,N-j,NumTrials)]),t,'support','positive','function','cdf');
        end
    end
end


end
