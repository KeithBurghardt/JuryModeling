function [InfluenceCDF] = GetInfluenceCDF(data,Tboundaries,Influencecdf,file)
    % Compute values using MLE within specified ranges
    Params = bInfluencefit(data,Tboundaries,Influencecdf,file);
    Binomialp_est = Params(1);
    MVMp_est      = Params(2);
    Tau_est       = Params(3);
    TauHung_est   = Params(4);
    Alpha_est     = Params(5);
    AlphaHung_est = Params(6);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%     Find ranges from file     %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [BinomialRange,MVMRange,TauRange,TauHungRange,AlphaRange,AlphaHungRange,N,~] = GetRanges(file);
    %{
    %string positions where "number of jurors" will be listed
    pos1=strfind(file,'N=');
    pos2=strfind(file,';Binomial_p='); 
    pos3=strfind(file,';MVM_p=');
    pos4=strfind(file,';Alpha=');

    pos5=strfind(file,';Alpha_hung=');
    pos6=strfind(file,';Tau=');
    pos7=strfind(file,';Tau_hung=');
    pos8=strfind(file,';NumRuns=');
    pos9=strfind(file,';NumTrials=');

    N = str2num(file(pos1+2:pos2-1));
    NumRuns = str2num(file(pos8+9:pos9-1));
    NumTrials = str2num(file(pos9+11:pos9+12));
    NumVals = NumRuns * NumTrials;

    %below order is the order in which the data is written
    BinomialRange=str2num(strrep(file(pos2+12:pos3-1),'_',':'));
    MVMRange=str2num(strrep(file(pos3+7:pos4-1),'_',':'));
    TauRange=str2num(strrep(file(pos6+5:pos7-1),'_',':'));
    TauHungRange=str2num(strrep(file(pos7+10:pos8-1),'_',':'));
    AlphaRange=str2num(strrep(file(pos4+7:pos5-1),'_',':'));
    AlphaHungRange=str2num(strrep(file(pos5+12:pos6-1),'_',':'));
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %N = numel(data(:,1)) - 1;
    % Influencecdf matrix position of parameters
    [~,Binomialp_pos] = ismember(Binomialp_est,BinomialRange);
    [~,MVMp_pos] = ismember(MVMp_est,MVMRange);
    [~,Tau_pos] = ismember(Tau_est,TauRange);     
    [~,TauHung_pos] = ismember(TauHung_est,TauHungRange);
    [~,Alpha_pos] = ismember(Alpha_est,AlphaRange);
    [~,AlphaHung_pos] = ismember(AlphaHung_est,AlphaHungRange);
    t_steps = 1:numel(Tboundaries);
    % CDF vs. vote corresponding to best fit parameters
    InfluenceCDF=reshape(Influencecdf(Binomialp_pos,MVMp_pos,Tau_pos,TauHung_pos,Alpha_pos,AlphaHung_pos,:,t_steps),N+1,numel(Tboundaries));
    
end
