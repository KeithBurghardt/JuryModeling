function [p, d, Dstar] = bTTpva(data,Tboundaries,TTcdf,LambdaRangeGuilty,LambdaRangeInnocent,varargin)

% bTTpva calculates the p-value for the given power-law fit to some data.
% Source: http://tuvalu.santafe.edu/~aaronc/powerlaws/bins/
%
%   When using binned data, the data vector 'h' is assumed to 
%   contain histogram counts between bin edges 'boundaries'.
%   Usage: a = plvar([900 90 9], [1 10 100 1000])
%   Note that while the above example uses logarithmic binning 
%   (powers of 10), any other binning scheme can be used in 
%   practice. 
%----------
% Options:
%----------  
% 1. a = bTTpva(h, boundaries, n, 'reps', 10000);
%    The default number of repetitions of fitting procedure is 
%    1000. This number can be changed using the reps option as shown
%    above
%
% 2. a = bTTpva(h, boundaries, n, 'silent');
%    This option can be used to silence the textual output on 
%    screen.
%
% See also BPLFIT, BPLVAR, BPLPVA, BEXPPVA
%
% Version 1.0 (2015) adapted from BPLPVA code:
% Copyright (C) 2012 Yogesh Virkar (University of Colorado, Boulder)/2015
% Keith Burghardt (University of Maryland, College Park)
% Distributed under GNU GPL v3.0
% http://www.gnu.org/copyleft/gpl.html
% bTTpva comes with ABSOLUTELY NO WARRANTY



% KS-test for 2 dimensions based on Fasano and Franceschini (1987)
% adapted from Press et al. "Numerical Recipes in C:
%                            The Art of Scientific Computing Second
%                            Edition", Cambridge Press, p.649

limit = [];
bminb = [];
reps = 1000;
silent = 0;

% ---------------------------------------------------------------
% ---------------Parsing command-line arguments------------------
% ---------------------------------------------------------------
i=1;
while i<=length(varargin)
    argok = 1;
    if(ischar(varargin{i}))
        switch varargin{i}
            case 'reps', reps = varargin{i+1}; i=i+1;
            case 'silent', silent = 1;
            otherwise, argok=0;    
        end
    end
    
    if ~argok,
        disp(['(bTTpva) Ignoring invalid argument #' num2str(i+2)]); 
    end
    i=i+1;
end


% ---------------------------------------------------------------
% ------------------------Checking input-------------------------
% ---------------------------------------------------------------

% 1. h must have integer counts.
if isequal(fix(reshape(data,numel(data),1)),reshape(data,numel(data),1))==0
    fprintf('(TTpva) Error: Vector h should be an integer vector.\n');
    return;
end

% 2. h must be non-negative
if ~isempty(find(reshape(data,numel(data),1)<0, 1))
    fprintf('(bTTpva) Error: Vector h should be non-negative.\n');
    return;
end

% 3. boundaries must have number of elements as one more than 
%    the number in h
if numel(Tboundaries)~=(numel(data(1,:))+1)
    fprintf('(bTTpva) Error: Incorrect number of elements in either boundaries or h.\n');
    return;
end

% 4. Need atleast 2 bins to work with.
if numel(data(1,:))<2
    fprintf('(bTTpva) Error: I need atleast 2 bins to make this work.\n');
    return;
end

% 5. Checking number of repititons
if ~isempty(reps) && (~isscalar(reps) || reps<2),
	fprintf('(bTTpva) Error: ''reps'' argument must be a positive value > 1; using default.\n');
    reps = 1000;
end;

% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------


d = zeros(reps,1);


if ~silent
    fprintf('TT consensus distribution, parameter uncertainty calculation\n');
    fprintf('    Copyright 2012 Yogesh Virkar/2015 Keith Burghardt\n');
    fprintf('    Warning: This can be a slow calculation; please be patient.\n');
    fprintf('    reps=%i\n', length(d));
end

% ---------------------------------------------------------------
%---------------Compute the empirical distance D*------------------
% ---------------------------------------------------------------




% Compute distance using KS statistic
TTCDF = GetTTCDF(data,Tboundaries,TTcdf,LambdaRangeGuilty,LambdaRangeInnocent);
Dstar=ks2d1sHist(data,TTCDF);
N = length(data(:,1)) - 1;
% ---------------------------------------------------------------
% Compute the distribution of gofs using parametric bootstrap
% ---------------------------------------------------------------

l = Tboundaries(1:end-1);
n = sum(reshape(data,numel(data),1));

parfor i=1:reps
    
    
    % create random numbers from TT
    %find the bootstrapped votes from the vote CDF
    z = TTCDF(:,end);
    temp = cumsum(z(end:-1:1));
    cx = 1 - temp(end:-1:1);
    % force cx(1) to be 0 (may be off by 10^-4)
    cx(1) = 0;
    v = Emprnd(cx,1:N+1,n);
    %NOTE: dimensions are incorrect
    hv = histc(v, 1:N+1);


    zq = zeros(N+1,length(Tboundaries) -1);
    for j = 1:N+1
        % bootstrapped histrogram of times for each vote
        % hv(j) is the number of votes (fraction of N total datapoints)
        % NOTE:TTCDF(j,end) < 1, thus we need to normalize
        newCDF = TTCDF(j,:)./TTCDF(j,end);
	newCDF(newCDF > 1) = 1;
        newCDF(:,1) = 0;
 	newCDF(newCDF < 0) = 0;

        x=Emprnd(newCDF,Tboundaries,hv(j));
        zq(j,:) = histc(x,Tboundaries(1:end-1));
    end
    bootCDF = GetTTCDF(zq,Tboundaries,TTcdf,LambdaRangeGuilty,LambdaRangeInnocent);
    dat = ks2d1sHist(zq,bootCDF);
  
    
    if ~silent
        fprintf('itr=%i\n', i);
    end
    d(i) = dat;
    
end
p = sum(d>=Dstar)./reps;


end
