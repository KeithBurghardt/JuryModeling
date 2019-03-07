function [CDF] = FindCDF_Influence(h,Influencecdf,Binomialp,MVMp,Tau,TauHung,Alpha,AlphaHung,file)

	[BinomialRange,MVMRange,TauRange,TauHungRange,AlphaRange,AlphaHungRange,N,NumVals] = GetRanges(file)
        [~,Binomialp_pos] = ismember(Binomialp,BinomialRange)
        [~,MVMp_pos] = ismember(MVMp,MVMRange)
        [~,Tau_pos] = ismember(Tau,TauRange)
        [~,TauHung_pos] = ismember(TauHung,TauHungRange)
        [~,Alpha_pos] = ismember(Alpha,AlphaRange)
        [~,AlphaHung_pos] = ismember(AlphaHung,AlphaHungRange)
        t_steps = 1:numel(h(1,:))+1;
        %make sure dimensions are correct
        N = numel(h(:,1))-1;
        fval = 0;
        CDF = zeros(N+1,numel(t_steps));
        for v = 1:N+1;
            % create binned PDF with the same dimensions as data
            temp = reshape(Influencecdf(Binomialp_pos,MVMp_pos,Tau_pos,TauHung_pos,Alpha_pos,AlphaHung_pos,v,t_steps),1,numel(t_steps));
            CDF(v,:) = temp;
        end

    end

