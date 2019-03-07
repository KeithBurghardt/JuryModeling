function [TTCDF] = GetTTCDF(data,Tboundaries,TTcdf,LambdaRangeGuilty,LambdaRangeInnocent)
% Compute values using MLE within specified ranges
    Params = bTTfit(data,Tboundaries,TTcdf,LambdaRangeGuilty,LambdaRangeInnocent);
    Lambda_Guilty = Params(1);
    Lambda_Innocent = Params(2);
    N = numel(data(:,1)) - 1;
    % TTcdf matrix position of parameters
    [~,lambda_guilty_pos]   = ismember(Lambda_Guilty,LambdaRangeGuilty);
    [~,lambda_innocent_pos] = ismember(Lambda_Innocent,LambdaRangeInnocent);
    t_steps = 1:numel(Tboundaries);
    % CDF vs. vote corresponding to best fit parameters
    TTCDF=reshape(TTcdf(lambda_guilty_pos,lambda_innocent_pos,:,t_steps),N+1,numel(Tboundaries));
    
    %TTCDF = TTcdf_paramspace(Lambda_Guilty,Lambda_Innocent,data, Tboundaries);
    %TTCDF = reshape(TTCDF,N+1,numel(Tboundaries));
end
