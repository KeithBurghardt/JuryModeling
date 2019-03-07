function [Sigma, eof] = bInfluencevar(data,Tboundaries,Influencecdf,file, varargin)

% bInfluencevar estimates the uncertainty in the estimated parameters.
%
%   When using binned data, the data vector 'h' is assumed to 
%   contain histogram counts between bin edges 'Tboundaries'.
%   Usage: a = bgevlvar([900 90 9], [1 10 100 1000])
%   Note that while the above example uses logarithmic binning 
%   (powers of 10), any other binning scheme can be used in 
%   practice. 
%
%----------
% Options:
%----------
% 1. a = bInfluencevar(h, Tboundaries, 'range', 1.5:0.01:3.5);
%    The 'range' option can be specified to restrict search for
%    T parameter. In above example, bgevlvar gives the best 
%    looking T in the specified range. By default bgevlvar uses 
%    matlab's fminsearch function which in turn uses the    
%    Nelder-Mead simplex search algorithm. Refer following url:
%    (http://www.mathworks.com/help/techdoc/math/bsotu2d.html#bsgpq6p-11)
%
% 2. a = bInfluencevar(h, Tboundaries, 'limit', 100);
%    The 'limit' option lets you limit the search for bmin. 
%    Values in Tboundaries above this limit are not considered as
%    candidate bmin values. 
%
% 3. a = bInfluencevar(h, Tboundaries, 'bmin', 100);
%    The 'bmin' option lets you fix a value such as 100 for bmin. 
%    Note that this value should be one of the values in the 
%    Tboundaries array. In the above example, 100 cannot be the 
%    last bin boundary. Also, it is advisable to give the fitting
%    procedure atleast two bins to work with. 
%
%    With options 2.and 3., if 'limit' or 'bmin' value is not one 
%    of the elements in Tboundaries, bInfluencevar chooses the bin boundary
%    which is closest to the specified value and less than that 
%    value.
%  
% 4. a = bInfluencevar(h, Tboundaries, 'reps', 10000);
%    The default number of repetitions of fitting procedure is 
%    1000. This number can be changed using the reps option as shown
%    above
%
% 5. a = bInfluencevar(h, Tboundaries, 'silent');
%    This option can be used to silence the textual output on 
%    screen.
%
% See also bInfluenceFIT, bInfluencePVA
%
% Version 1.0 (2012)
% Copyright (C) 2012 Yogesh Virkar (University of Colorado, Boulder)/2016
% Keith Burghardt (University of Maryland, College Park)
% Distributed under GNU GPL v3.0
% http://www.gnu.org/copyleft/gpl.html
% bInfluencevar comes with ABSOLUTELY NO WARRANTY

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
        disp(['(bInfluencevar) Ignoring invalid argument #' num2str(i+2)]); 
    end
    i=i+1;
end

% ---------------------------------------------------------------
% ------------------------Checking input-------------------------
% ---------------------------------------------------------------


% 1. Tboundaries must have number of elements as one more than 
%    the number in h

if numel(Tboundaries)~=(numel(data(1,:))+1)
    fprintf('(bInfluencevar) Error: Incorrect number of elements in either Tboundaries or h.\n');
    return;
end

% 2. Need atleast 2 bins to work with.
if numel(data(1,:))<2
    fprintf('(bInfluencevar) Error: I need atleast 2 bins to make this work.\n');
    return;
end

% 3. Checking range vector
if ~isempty(rngal) && (~isvector(rngal) || min(rngal)<1)
    fprintf('(bInfluencevar) Error: ''range'' argument must contain a valid vector; using default.\n');
    rngal = 1.5:0.01:3.5;
end

% 4. Checking limit option
if ~isempty(limit) && (~isscalar(limit) || limit<min(Tboundaries))
    fprintf('(bInfluencevar) Error: ''limit'' argument must be a positive value >= Tboundaries(1); using default.\n');
    limit = Tboundaries(end-2);
end

% 5. Checking bmin option
if ~isempty(bminb) && (~isscalar(bminb) || bminb>=Tboundaries(end-1))
    fprintf('(bInfluencevar) Error: ''bmin'' argument must be a positive value < Tboundaries(end-1); using default.\n');
    bminb = Tboundaries(1);
end

% 6. Checking number of repititons
if ~isempty(reps) && (~isscalar(reps) || reps<2),
	fprintf('(bInfluencevar) Error: ''reps'' argument must be a positive value > 1; using default.\n');
    reps = 1000;
end;

% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------


N = numel(data(:,1)) - 1;
n = sum(reshape(data,numel(data),1));
eof = zeros(reps,6);

if ~silent
    fprintf('Influence consensus distribution, parameter uncertainty calculation\n');
    fprintf('    Copyright 2012 Yogesh Virkar/2016 Keith Burghardt\n');
    fprintf('    Warning: This can be a slow calculation; please be patient.\n');
    fprintf('    reps=%i\n', length(eof(:,1)));
end

parfor i=1:reps
    
    % Bootstrap resample of binned data
    l = Tboundaries(1:end-1);
    u = Tboundaries(2:end);
    temp = (l+u)./2;
    Vihist = sum(data,2);
    temp2 = cumsum(Vihist(end:-1:1));
    cx = 1 - temp2(end:-1:1)./n;

    Nviboot = histc(Emprnd(cx,1:N+1,n),1:N+1);

    bootdata = zeros(size(data));
    
    % for vote vi in data
    for vi = 1:N + 1
        temp2=[];
        % for all elements in data with vote vi
        for y=1:numel(data(vi,:))
            % repeat all data data(vi,y) times for data value = temp(y)
            temp2 = [temp2;repmat(temp(y),round(data(vi,y)),1)];
        end
        
        % total number of data points
        Ni = round(sum(data(vi,:)));
        %randomly permute values
        temp2 = temp2(randperm(Ni));
        % bootstraped points: Niboot values with position 1:Ni
        xboot = temp2(ceil(Ni*rand(Nviboot(vi),1)));
        %make new histogram
        [h2, ~] = histc(xboot, Tboundaries);
        bootdata(vi,:) = h2(1:end - 1);
    end

    % estimate T using specified range
    Params = bInfluencefit(bootdata,Tboundaries,Influencecdf,file);
    Binomialp_est = Params(1);
    MVMp_est      = Params(2);
    Tau_est       = Params(3);
    TauHung_est   = Params(4);
    Alpha_est     = Params(5);
    AlphaHung_est = Params(6);
    % ni = sum(reshape(bootdata,numel(bootdata),1));
    % Fraction of trials that vote guilty
    % Store all parameter values.
    eof(i,:) = Params;
    if ~silent
        fprintf('itr=%i, Binomial p=%f, MVM p=%f, Tau=%f, Tau hung=%f, Alpha=%f, Alpha hung=%f\n', i,Binomialp_est,MVMp_est,Tau_est,TauHung_est,Alpha_est,AlphaHung_est);
    end
    
end

% Calculate the uncertainty
sigma_Binomialp = std(eof(:,1));
sigma_MVMp      = std(eof(:,2));
sigma_Tau       = std(eof(:,3));
sigma_TauHung   = std(eof(:,4));
sigma_Alpha     = std(eof(:,5));
sigma_AlphaHung = std(eof(:,6));
Sigma = [sigma_Binomialp sigma_MVMp sigma_Tau sigma_TauHung sigma_Alpha sigma_AlphaHung];

end
