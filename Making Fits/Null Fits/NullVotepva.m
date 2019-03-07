function [p, d, Dstar,Params] = NullVotepva(h, varargin)

% BEXPPVA calculates the p-value for the given exponential fit to some data.
% Source: http://tuvalu.santafe.edu/~aaronc/powerlaws/bins/
%
%   When using binned data, the data vector 'h' is assumed to 
%   contain histogram counts between bin edges 'boundaries'.
%   Usage: a = plvar([900 90 9], [1 10 100 1000])
%   Note that while the above example uses logarithmic binning 
%   (powers of 10), any other binning scheme can be used in 
%   practice. 
%
%----------
% Options:
%----------  
% 1. a = bexppva(h, boundaries, bmin, 'reps', 10000);
%    The default number of repetitions of fitting procedure is 
%    1000. This number can be changed using the reps option as shown
%    above
%
% 2. a = bexppva(h, boundaries, bmin, 'silent');
%    This option can be used to silence the textual output on 
%    screen.
%
% See also BPLFIT, BPLVAR, BPLPVA
%
% Version 1.0 (2015) adapted from BPLPVA code:
% Copyright (C) 2012 Yogesh Virkar (University of Colorado, Boulder)
% Distributed under GNU GPL v3.0
% http://www.gnu.org/copyleft/gpl.html
% BPLPVA/BEXPPVA comes with ABSOLUTELY NO WARRANTY


rngal = [];
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
            case 'range', rngal = varargin{i+1}; i=i+1;
            case 'limit', limit = varargin{i+1}; i=i+1;
            case 'bmin', bminb = varargin{i+1}; i=i+1;
            case 'reps', reps = varargin{i+1}; i=i+1;
            case 'silent', silent = 1;
            otherwise, argok=0;    
        end
    end
    
    if ~argok,
        disp(['(BEXPPVA) Ignoring invalid argument #' num2str(i+2)]); 
    end
    i=i+1;
end


% ---------------------------------------------------------------
% ------------------------Checking input-------------------------
% ---------------------------------------------------------------

% 1. h must have integer counts.
if isequal(fix(h),h)==0
    fprintf('(BEXPPVA) Error: Vector h should be an integer vector.\n');
    return;
end

% 2. h must be non-negative
if ~isempty(find(h<0, 1))
    fprintf('(BEXPPVA) Error: Vector h should be non-negative.\n');
    return;
end


% 4. Need atleast 2 bins to work with.
if numel(h)<2
    fprintf('(BEXPPVA) Error: I need atleast 2 bins to make this work.\n');
    return;
end

% 5. Checking range vector
if ~isempty(rngal) && (~isvector(rngal) || min(rngal)<1)
    fprintf('(BEXPPVA) Error: ''range'' argument must contain a valid vector; using default.\n');
    rngal = 1.5:0.01:3.5;
end


% 8. Checking number of repititons
if ~isempty(reps) && (~isscalar(reps) || reps<2),
	fprintf('(BEXPPVA) Error: ''reps'' argument must be a positive value > 1; using default.\n');
    reps = 1000;
end;

% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------

% Reshape the input vectors
h = reshape(h, 1,numel(h));
%Number of jurors
n = numel(h) - 1;
votes = 0:n;

N = sum(h);
d = zeros(reps,1);


if ~silent
    fprintf('Exponential distribution, parameter uncertainty calculation\n');
    fprintf('    Copyright 2012 Yogesh Virkar/2015 Keith Burghardt\n');
    fprintf('    Warning: This can be a slow calculation; please be patient.\n');
    fprintf('    reps=%i\n', length(d));
end

% ---------------------------------------------------------------
%---------------Compute the empirical distance D*------------------
% ---------------------------------------------------------------
% find the cdf for the data

NumInnocent = sum((0:floor(n/2)).*h(1:floor(n/2)+1));
MeanInnocent = NumInnocent/sum(h(1:floor(n/2)+1));
PInnocent = MeanInnocent/n;

NumGuilty = sum((floor(n/2)+1:n).*h(floor(n/2)+2:n+1));
MeanGuilty = NumGuilty./sum(h(floor(n/2)+2:n+1));
PGuilty = MeanGuilty/n;
Params = [PInnocent PGuilty];
pf = pdf('Binomial',votes,n,PInnocent).*sum(h(1:floor(n/2)+1))/N + pdf('Binomial',votes,n,PGuilty).*sum(h(floor(n/2)+2:n+1))/N;
temp = cumsum(pf(end:-1:1));
cf = 1 - temp(end:-1:1);
BestFitModel = cf;

% Compute lambda for exponential fit%disp(h);
%disp(boundaries);
% Compute distance using KS statistic
temp = cumsum(h(end:-1:1));
cx = 1 - temp(end:-1:1)./N;


%figure;
%plot(1:length(cx),cx);hold on;plot(1:length(cf),cf,'r');

Dstar = max(abs(cf-cx));

% ---------------------------------------------------------------
% Compute the distribution of gofs using semiparametric bootstrap
% ---------------------------------------------------------------

%number of repetitions, subtracting when the histogram has 1 bin
numreps = reps;

parfor i=1:reps
    %semi-parametric bootstrap of data
    temp = votes;
    x = Emprnd(BestFitModel,votes,N)+0.5;
    h2 = histc(x, [votes n+1]);
    h2 = h2(votes + 1);
    % Compute distance using KS statistic
    temp = cumsum(h2(end:-1:1));
    cx = 1 - temp(end:-1:1)./N;

    NumInnocent = sum((0:floor(n/2)).*h2(1:floor(n/2)+1));
    MeanInnocent = NumInnocent/sum(h2(1:floor(n/2)+1));
    PInnocent = MeanInnocent/n;
    
    NumGuilty = sum((floor(n/2)+1:n).*h2(floor(n/2)+2:n+1));
    MeanGuilty = NumGuilty./sum(h2(floor(n/2)+2:n+1));
    PGuilty = MeanGuilty/n;
    pf = pdf('Binomial',votes,n,PInnocent).*sum(h2(1:floor(n/2)+1))/N + pdf('Binomial',votes,n,PGuilty).*sum(h2(floor(n/2)+2:n+1))/N;
    temp = cumsum(pf(end:-1:1));
    cf = 1 - temp(end:-1:1);

    if ~silent
        fprintf('itr=%i\n', i);
    end
    d(i) = max(abs(cf-cx));

end

p = sum(d>=Dstar)./numreps;

end
