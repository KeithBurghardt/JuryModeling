function [lambda_est, L] = bexpnfit(h, boundaries, varargin)

% BEXPFIT fits a exponential model to the data
% Source: http://tuvalu.santafe.edu/~aaronc/powerlaws/bins/
%
%----------
% Options:
%----------
% 1. l = bexpnfit(h, boundaries, 'range', 0.1:0.01:1)
%    The 'range' option can be specified to restrict search for
%    lambda parameter. In above example, bexpnfit gives the best 
%    looking alpha in the specified range. By default bexpnfit uses 
%    matlab's fminsearch function which in turn uses the    
%    Nelder-Mead simplex search algorithm. Refer following url:
%    (http://www.mathworks.com/help/techdoc/math/bsotu2d.html#bsgpq6p-11)
%
% 2. l = bexpnfit(h, boundaries, 'bmin', 100)
%    The 'bmin' option lets you fix a value such as 100 for bmin. 
%    Note that this value should be one of the values in the 
%    boundaries array. In the above example, 100 cannot be the 
%    last bin boundary. Also, it is advisable to give the fitting
%    procedure atleast two bins to work with. 
%    Note that if 'bmin' value is not one of the elements in 
%    boundaries, bexpnfit chooses the bin boundary which is closest 
%    to the specified value and less than that value. Also, the 
%    Default bmin value is the first bin boundary i.e. 
%    boundaries(1).
%
% Version 1.0 (2012)
% Copyright (C) 2012 Yogesh Virkar (University of Colorado, Boulder)
% Distributed under GNU GPL v3.0
% http://www.gnu.org/copyleft/gpl.html
% BEXPFIT comes with ABSOLUTELY NO WARRANTY

rngal = [];
bminb = [];

% ---------------------------------------------------------------
% ---------------Parsing command-line arguments------------------
% ---------------------------------------------------------------
i=1;
while i<=length(varargin)
    argok = 1;
    if(ischar(varargin{i}))
        switch varargin{i}
            case 'range', rngal = varargin{i+1}; i=i+1;
            case 'bmin', bminb = varargin{i+1}; i=i+1;
            otherwise, argok=0;    
        end
    end
    if ~argok,
        disp(['(BEXPFIT) Ignoring invalid argument #' num2str(i+2)]); 
    end
    i=i+1;
end

% ---------------------------------------------------------------
% ------------------------Checking input-------------------------
% ---------------------------------------------------------------

% 1. h must have integer counts.
if isequal(fix(h),h)==0
    fprintf('(BEXPFIT) Error: Vector h should be an integer vector.\n');
    return;
end

% 2. h must be non-negative
if ~isempty(find(h<0, 1))
    fprintf('(BEXPFIT) Error: Vector h should be non-negative.\n');
    return;
end

% 3. boundaries must have number of elements as one more than 
%    the number in h
if numel(boundaries)~=(numel(h)+1)
    fprintf('(BEXPFIT) Error: Incorrect number of elements in either boundaries or h.\n');
    return;
end

% 4. Need atleast 2 bins to work with.
if numel(h)<2
    fprintf('(BEXPFIT) Error: I need atleast 2 bins to make this work.\n');
    return;
end

% 5. Checking range vector
if ~isempty(rngal) && (~isvector(rngal) || min(rngal)<=0)
    fprintf('(BEXPFIT) Error: ''range'' argument must contain a valid vector; using default.\n');
    rngal = 0.1:0.01:1;
end

% 6. Checking bmin option
if ~isempty(bminb) && (~isscalar(bminb) || bminb>=boundaries(end-1))
    fprintf('(BEXPFIT) Error: ''bmin'' argument must be a positive value < boundaries(end-1); using default.\n');
    bminb = [];
end

% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------

% Reshape the input vectors
h = reshape(h, numel(h), 1);
boundaries = reshape(boundaries, numel(boundaries), 1);

% Default bmin used is boundaries(1)
if isempty(bminb)
    bminb = boundaries(1);
end

ind = find(boundaries<=bminb, 1, 'last');
bminb = boundaries(ind);
h2 = h(ind:end);
boundaries2  = boundaries(ind:end);

l = boundaries2(1:end-1);
u = boundaries2(2:end);

    function [fval] = exponential_mle(lambda)
        temp = exp(-lambda.*l) - exp(-lambda.*u);
        %temp =  exp(-lambda).*(exp(l) - exp(u));
        temp = temp + 0.000001;
        fval = -sum(h2.*log(temp) + bminb.*(lambda.*h2)); 
    end

if ~isempty(rngal)
    rngal = reshape(rngal, 1, numel(rngal));
    LAMBDA = repmat(rngal, numel(l), 1);
    disp(LAMBDA);
    L = repmat(l, 1, numel(rngal));
    disp(L);
    U = repmat(u, 1, numel(rngal));
    H = repmat(h2, 1, numel(rngal));
    
    temp = exp(-LAMBDA.*L) - exp(-LAMBDA.*U);
    temp = temp + 0.000001;
    temp2  = H.*log(temp) + bminb.*(H.*LAMBDA);
    fval  = -sum(temp2, 1);
    [~,I] = min(fval);
    lambda_est = rngal(I);
    plot(fval);
    
else
    %options = optimset('MaxIter', 100);
    %lambda_est = fminsearch(@exponential_mle, 0.01, options);
    lambda_est = fminsearchbnd(@exponential_mle, 0.01,0,10);
end

% Following condition typically occurs when EXP is a bad fit
if(lambda_est>20)
    lambda_est = 20;
%elseif(lambda_est < 0)
%    lambda_est = 0.00001;
end

L = -sum(h2.*log(exp(-lambda_est.*l) - exp(-lambda_est.*u)) + (bminb*lambda_est).*h2);

end