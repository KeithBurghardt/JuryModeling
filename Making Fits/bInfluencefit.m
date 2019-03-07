function [param_est, L] = bInfluencefit(data,Tboundaries,Influencecdf,file)

% BInfluenceFIT fits a Influence model to the data from a sampled distribution
%
% Copyright (C) 2015 Keith Burghardt (University of Maryland, College Park) based on code
% by  Yogesh Virkar (University of Colorado, Boulder)
% Distributed under GNU GPL v3.0
% hInfluencep://www.gnu.org/copyleft/gpl.html
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
%  FiInfluenceing model to data: 
%       Maximize the Likelihood function/Log-likelihood function
%      - Likelihood function     L = Prod_i p(ti,vi|params)
%      - Log-likelihood function l = Sum_i log(p(ti,vi|params))
%       
%
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%     Find ranges from file     %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%t = reshape(Tboundaries,numel(Tboundaries),1);
[BinomialRange,MVMRange,TauRange,TauHungRange,AlphaRange,AlphaHungRange,N,~] = GetRanges(file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape the input vectors
Tboundaries = reshape(Tboundaries, numel(Tboundaries), 1);
l = Tboundaries(1:end-1);
u = Tboundaries(2:end);
%N = length(data(:,1)) - 1;

%n = sum(reshape(data,numel(data),1));

    function [fval] = Influence_mle(Binomialp,MVMp,Tau,TauHung,Alpha,AlphaHung)
        % CDF (T,t)
        %FullModelCDF = Influencecdf_paramspace(Tau_Guilty,Tau_Innocent,data,Tboundaries,N);
        [~,Binomialp_pos] = ismember(Binomialp,BinomialRange);
        [~,MVMp_pos] = ismember(MVMp,MVMRange);
        [~,Tau_pos] = ismember(Tau,TauRange);
        [~,TauHung_pos] = ismember(TauHung,TauHungRange);
        [~,Alpha_pos] = ismember(Alpha,AlphaRange);
        [~,AlphaHung_pos] = ismember(AlphaHung,AlphaHungRange);
        l_steps = 1:numel(data(1,:));
        u_steps = 2:numel(data(1,:))+1;
        %make sure dimensions are correct
        
        fval = 0;
        %temp = Influencecdf(lambda_guilty_pos,lambda_innocent_pos,:,u_steps) - Influencecdf(lambda_guilty_pos,lambda_innocent_pos,:,l_steps);
        for v = 1:N+1;
            % create binned PDF with the same dimensions as data
            
            temp = reshape(Influencecdf(Binomialp_pos,MVMp_pos,Tau_pos,TauHung_pos,Alpha_pos,AlphaHung_pos,v,u_steps)...
                - Influencecdf(Binomialp_pos,MVMp_pos,Tau_pos,TauHung_pos,Alpha_pos,AlphaHung_pos,v,l_steps),1,length(l_steps));
            
            temp = temp + 10^-11;%0.000001;
            fval = fval - sum(data(v,:).*log(temp));% CHECK
        end

    end
    function [L] = Influence_getL(b_index,mvm_index,tau_index,tau_hung_index,a_index,a_hung_index)

        l_steps = 1:numel(Tboundaries) - 1;
        u_steps = 2:numel(Tboundaries);
	L = zeros(N+1,sum(sum(data)));
	count = 1;
        for v = 1:N+1;
            p = reshape(Influencecdf(b_index,mvm_index,tau_index,tau_hung_index,a_index,a_hung_index,v,u_steps) ...
		- Influencecdf(b_index,mvm_index,tau_index,tau_hung_index,a_index,a_hung_index,v,l_steps),1,numel(l_steps));

            temp = (l+u)./2;
            temp2=[];
            for y=1:numel(data(v,:))
                temp2 = [temp2;repmat(temp(y),data(v,y),1)];
            end
	    N2 = numel(temp2);
            temp2 = temp2(randperm(N2));
            [~,whichbin] = histc(temp2, Tboundaries);

            % Calculate log-likelihood at each data point.
            Lv = log(p(whichbin));

            NEGINF = -100000000;
            Lv(Lv==-Inf) = NEGINF;
	    L(count:count+numel(Lv)-1) = Lv;
	    count = count + numel(Lv);
	end
    end

fval = zeros(numel(BinomialRange),numel(MVMRange),numel(TauRange),numel(TauHungRange),numel(AlphaRange),numel(AlphaHungRange));
for b = 1:numel(BinomialRange)
    for mvm = 1:numel(MVMRange)
        for tau = 1:numel(TauRange)
            for tau_hung = 1:numel(TauHungRange)
                for a = 1:numel(AlphaRange)
                    for a_hung = 1:numel(AlphaHungRange)
                        Binomialp = BinomialRange(b);
                        MVMp      = MVMRange(mvm);
                        Tau       = TauRange(tau);
                        TauHung   = TauHungRange(tau_hung);
                        Alpha     = AlphaRange(a);
                        AlphaHung = AlphaHungRange(a_hung);
                        fval(b,mvm,tau,tau_hung,a,a_hung)=Influence_mle(Binomialp,MVMp,Tau,TauHung,Alpha,AlphaHung);
                    end
                end
            end
        end
    end
end
[~,I] = min(reshape(fval,numel(fval),1));
[b_index,mvm_index,tau_index,tau_hung_index,a_index,a_hung_index] = ind2sub(size(fval),I);
% find positions
Binomialp_est = BinomialRange(b_index);
MVMp_est      = MVMRange(mvm_index);
Tau_est       = TauRange(tau_index);
TauHung_est   = TauHungRange(tau_hung_index);
Alpha_est     = AlphaRange(a_index);
AlphaHung_est = AlphaHungRange(a_hung_index);

param_est = [Binomialp_est MVMp_est Tau_est TauHung_est Alpha_est AlphaHung_est];
L = Influence_getL(b_index,mvm_index,tau_index,tau_hung_index,a_index,a_hung_index);

end
