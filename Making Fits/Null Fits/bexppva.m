function [p, d, Dstar] = bexppva(h, boundaries, bmin, varargin)

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

% 3. boundaries must have number of elements as one more than 
%    the number in h
if numel(boundaries)~=(numel(h)+1)
    fprintf('(BEXPPVA) Error: Incorrect number of elements in either boundaries or h.\n');
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

% 6. Checking limit option
if ~isempty(limit) && (~isscalar(limit) || limit<min(boundaries))
    fprintf('(BEXPPVA) Error: ''limit'' argument must be a positive value >= boundaries(1); using default.\n');
    limit = boundaries(end-2);
end

% 7. Checking bmin option
if ~isempty(bminb) && (~isscalar(bminb) || bminb>=boundaries(end-1))
    fprintf('(BEXPPVA) Error: ''bmin'' argument must be a positive value < boundaries(end-1); using default.\n');
    bminb = boundaries(1);
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
h = reshape(h, numel(h), 1);
boundaries = reshape(boundaries, numel(boundaries), 1);

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
% Data above bmin
ind = find(boundaries>=bmin, 1);
z = h(ind:end);     nz = sum(z);
b = boundaries(ind:end);
l = b(1:end-1);
u = b(2:end);

% Data below bmin
y = h(1:ind-1);     %ny = sum(y);
by = boundaries(1:ind);
ly = by(1:end-1);
uy = by(2:end);

% Compute alpha using numerical maximization
%hnd = @(alpha) -sum( z.*( log((l).^(1-alpha) - (u).^(1-alpha)) + (alpha-1)*log(bmin) ) );
%alpha = fminsearch(hnd, 1);
% Compute lambda for exponential fit
lambda=bexpnfit(h,boundaries,'bmin',bmin);
%disp(h);
%disp(boundaries);
% Compute distance using KS statistic
temp = cumsum(z(end:-1:1));
cx = 1 - temp(end:-1:1)./nz;
%cf = 1 - (l./bmin).^(1-alpha); 
cf = 1 - exp(-lambda*(l-bmin));

%figure;
%plot(1:length(cx),cx);hold on;plot(1:length(cf),cf,'r');

Dstar = max(abs(cf-cx));    

% ---------------------------------------------------------------
% Compute the distribution of gofs using semiparametric bootstrap
% ---------------------------------------------------------------

% Probability of choosing value above bmin 
pz = nz/N;
%number of repetitions, subtracting when the histogram has 1 bin
numreps = reps;

parfor i=1:reps
    %semi-parametric bootstrap of data
    n1 = sum(rand(1,N)>pz);
    temp = (ly+uy)./2;
    temp2=[];
    for t=1:numel(y)
        temp2 = [temp2;repmat(temp(t),y(t),1)];
    end
    temp2 = temp2(randperm(numel(temp2)));
    x1 = temp2(ceil(numel(temp2)*rand(n1,1)));
    n2 = N-n1;
    %x2 = bmin.*(1-rand(n2,1)).^(-1/(alpha-1));
    %n2 random numbers for an exponential PDF
    %NOTE: minus sign
    x2 = bmin + exprnd(1/lambda,n2,1);%-bmin.*log(1-rand(n2,1))./(lambda);
    x = [x1;x2];
    h2 = histc(x, boundaries);
    h2(end)=[];
    ind = find(h2(end:-1:1)~=0,1,'first')-1;
    if(ind==1)
        h2(end)= [];
    else
        if(ind>1)
            ind2=ind-1;
            h2(end-ind2:end) = [];
        end    
    end
        
    boundaries2 = boundaries(1:end-ind);
    
    % Need a minimum of 2 bins.
    bmins = boundaries2(1:end-2);
    
    if ~isempty(bminb)
        bmins = bmins(find(bmins<=bminb, 1, 'last'));
    end
    if ~isempty(limit)
        bmins(bmins>limit) = [];
    end

    dat = zeros(size(bmins));
    for xm=1:length(bmins)
        %if(lam >= -999 && lam ~=20)    
            bminq = bmins(xm);
    
            % Truncate the data below bmin
            indq = find(boundaries2>=bminq, 1);
            zq = h2(indq:end);
            nq = sum(zq);
            bq = boundaries2(indq:end); 
        
            % estimate alpha using specified range or using 
            % numerical maximization
            lq = bq(1:end-1);
            uq = bq(2:end);
            % estimate lambda using numerical maximization
            %h3 = histc(x, boundaries);
            %h3(end)=[];
            
            lam=bexpnfit(zq,bq,'bmin',bminq);
            %disp(lam);
            % compute KS statistic
            tempq = cumsum(zq(end:-1:1));
            cxq = 1 - tempq(end:-1:1)./nq;
            %cfq = 1 - (lq./bminq).^(1-al); 
            cfq = 1 - exp(-lam.*(lq-bminq));
            dat(xm) = max(abs(cfq-cxq));
        %else
        %    dat(xm) = 999999;
        %end
    end
    
    if ~silent
        fprintf('itr=%i\n', i);
    end
    if ~isempty(dat)
        d(i) = min(dat);
    else
        %this value is not counted, because histogram was too poor
        d(i) = 0;
        numreps = numreps - 1;
    end
end


p = sum(d>=Dstar)./numreps;


end