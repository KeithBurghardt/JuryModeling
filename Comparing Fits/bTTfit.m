function [param_est, L,sumL] = bTTfit(data,Tboundaries,TTcdf,LambdaRangeGuilty,LambdaRangeInnocent)

% BTTFIT fits a TT model to the data from a sampled distribution
%
% Copyright (C) 2015 Keith Burghardt (University of Maryland, College Park) based on code
% by  Yogesh Virkar (University of Colorado, Boulder)
% Distributed under GNU GPL v3.0
% http://www.gnu.org/copyleft/gpl.html
% BFULLMODELFIT comes with ABSOLUTELY NO WARRANTY
%
%   Data:  - Of the form data(i,:) = {ti:final decision is vi}
%           e.g.,     (vi = 0-12) [ 0 10 1023 203 100 10 2 0 ]
%                     (vi = 1-11) [ 10 120 123 23 1 0 0 0    ]
%                     (vi = 2-10) [ 4 11 103 1203 120 10 1 16]
%                     (vi = 3-9 ) [ 2 210 18 22 13 12 1 0    ]
%                       ...         ...
%
%   Model: - Params:
%
%               1. Lambda_Guilty:    the timescale to stop if voting guilty 
%                                       
%               2. Lambda_Innocent:  the timescale to stop if voting innocent          
%
%  Fitting model to data: 
%       Maximize the Likelihood function/Log-likelihood function
%      - Likelihood function     L = Prod_i p(ti,vi|params)
%      - Log-likelihood function l = Sum_i log(p(ti,vi|params))
%       
%
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------

% Reshape the input vectors
Tboundaries = reshape(Tboundaries, numel(Tboundaries), 1);
%l = Tboundaries(1:end-1);
%u = Tboundaries(2:end);
N = length(data(:,1)) - 1;

%n = sum(reshape(data,numel(data),1));

    function [fval] = TT_mle(Lambda_Guilty,Lambda_Innocent)
        % CDF (T,t)
        %FullModelCDF = TTcdf_paramspace(Tau_Guilty,Tau_Innocent,data,Tboundaries,N);
        [~,lambda_guilty_pos] = ismember(Lambda_Guilty,LambdaRangeGuilty);
        [~,lambda_innocent_pos] = ismember(Lambda_Innocent,LambdaRangeInnocent);
        % H = repmat(data,length(lambda),1);
        l_steps = 1:numel(data(1,:));
        u_steps = 2:numel(data(1,:))+1;
        %make sure dimensions are correct
        temp = 0; 
        fval = 0;
        %temp = TTcdf(lambda_guilty_pos,lambda_innocent_pos,:,u_steps) - TTcdf(lambda_guilty_pos,lambda_innocent_pos,:,l_steps);
        for v = 1:N+1;
            % create binned PDF with the same dimensions as data
            
            temp = reshape(TTcdf(lambda_guilty_pos,lambda_innocent_pos,v,u_steps) - TTcdf(lambda_guilty_pos,lambda_innocent_pos,v,l_steps),1,length(l_steps));
            
            temp = temp + 0.000001;
            fval = fval - sum(data(v,:).*log(temp));% CHECK
        end

    end
    
    function [L] =  TT_getL(Lambda_Guilty,Lambda_Innocent)
        [~,lambda_guilty_pos] = ismember(Lambda_Guilty,LambdaRangeGuilty);
        [~,lambda_innocent_pos] = ismember(Lambda_Innocent,LambdaRangeInnocent);
        
        min_pos = 1;%find(Xrange == l(1));
        l_steps = min_pos:min_pos + numel(data(1,:)) - 1;
        u_steps = min_pos + 1:min_pos + numel(data(1,:)); 
        p = TTcdf(lambda_guilty_pos,lambda_innocent_pos,:,u_steps) - TTcdf(lambda_guilty_pos,lambda_innocent_pos,:,l_steps);
        
        temp = (l_steps+u_steps)./2;
        temp2=[];
        %fix
        for i = 1:N+1
            for y=1:numel(data(i,:))
                temp2 = [temp2;repmat(temp(y),data(y),1)];
            end
            N = numel(temp2);    
            temp2 = temp2(randperm(N));
            [~,whichbin] = histc(temp2,min_pos:min_pos + numel(data(1,:)));%Tboundaries);
            % Calculate log-likelihood at each data point.
            L = log(p(whichbin));
        end
        

        NEGINF = -100000000;
        L(L==-Inf) = NEGINF;

    end
    
%param_est = fminsearchbnd(@(Tau) FullModel_mle(Tau(1),Tau(2)), [0.01 0.01],[0 0],[10 10]);
%fval = TT_mle(LambdaRangeGuilty,LambdaRangeInnocent);
fval = zeros(length(LambdaRangeGuilty),length(LambdaRangeInnocent));
%takes a LOOOOONG time. See if this can be sped up
counti = 0;
for i = LambdaRangeGuilty
	counti = counti + 1;
	countj = 0;
	for j = LambdaRangeInnocent
		countj = countj + 1;
		fval(counti,countj)=TT_mle(i,j);
	end
end
[sumL,I] = min(reshape(fval,numel(fval),1));
[lambda_guilty_index,lambda_innocent_index] = ind2sub(size(fval),I);
% find positions
Lambda_Guilty_est = LambdaRangeGuilty(lambda_guilty_index);%param_est(1);
Lambda_Innocent_est= LambdaRangeInnocent(lambda_innocent_index);%param_est(2);
param_est = [Lambda_Guilty_est Lambda_Innocent_est];

L = TT_getL(Lambda_Guilty_est,Lambda_Innocent_est);
%L=0;

end
