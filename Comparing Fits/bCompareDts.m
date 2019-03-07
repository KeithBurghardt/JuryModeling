function [R,p] = bCompareDts(h, Tboundaries,fileDt15,fileDt60,fileDt240)
    % What this function does:
    %   - Find estimates of values, error bars, & bmin (if exponential) 
    %   - Find p-value of each function
    %   - If p-value > 0.1, compare noramalized likelihood ratios
    %

    %Total amount of data
    n = sum(reshape(h,numel(h),1));

    % number of jurors
    N = length(h(:,1)) - 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   Majority Voter Model - Stubbornness (3)   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = 15;
    N = numel(h(:,1)) - 1;
    %if N == 6
    %file = '2StrainJuryDelib-3sec;Beta=1;N=6;Binomial_p=0.53;MVM_p=0.88_0.02_1.00;Alpha=0.0002_0.0002_0.002;Alpha_hung=0.00003_0.00003_0.0001;Tau=20_80_420;Tau_hung=0.05_0.05_0.15;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    %if N == 8
	%file = '2StrainJuryDelib-3sec;Beta=1;N=8;Binomial_p=0.53;MVM_p=0.6_0.035_0.8;Alpha=0.0002_0.0002_0.002;Alpha_hung=0.00002_0.00002_0.0001;Tau=40_60_400;Tau_hung=0.02_0.04_0.15;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    %elseif N == 12
	%file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0.53;MVM_p=0.65_0.04_0.97;Alpha=0.0001_0.0002_0.002;Alpha_hung=0.00002_0.00003_0.00015;Tau=20_60_400;Tau_hung=0.02_0.05_0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    %end

    InfluenceCDF = Influencecdf_paramspace(fileDt15,dt,Tboundaries);
    [~, lMVM_15] = bInfluencefit(h, Tboundaries,InfluenceCDF,fileDt15);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   Majority Voter Model - Stubbornness (60)   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = 60;
    N = numel(h(:,1)) - 1;
    %if N == 6
    %file = '2StrainJuryDelib-3sec;Beta=1;N=6;Binomial_p=0.53;MVM_p=0.88_0.02_1.00;Alpha=0.004_0.004_0.04;Alpha_hung=0.0006_0.0006_0.002;Tau=1_4_21;Tau_hung=0.05_0.05_0.15;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    %if N == 8
    %    file = '2StrainJuryDelib-3sec;Beta=1;N=8;Binomial_p=0.53;MVM_p=0.6_0.035_0.8;Alpha=0.004_0.004_0.04;Alpha_hung=0.0004_0.0004_0.002;Tau=2_3_20;Tau_hung=0.02_0.04_0.15;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    %elseif N == 12
    %	file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0.53;MVM_p=0.65_0.04_0.97;Alpha=0.002_0.004_0.04;Alpha_hung=0.0004_0.0006_0.003;Tau=1_3_20;Tau_hung=0.02_0.05_0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    %end


    InfluenceCDF = Influencecdf_paramspace(fileDt60,dt,Tboundaries);
    [~, lMVM_60] = bInfluencefit(h, Tboundaries,InfluenceCDF,fileDt60);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   Majority Voter Model - Stubbornness (300)   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = 240;
    N = numel(h(:,1)) - 1;
    %file = '2StrainJuryDelib-3sec;Beta=1;N=6;Binomial_p=0.53;MVM_p=0.88_0.02_1.00;Alpha=0.02_0.02_0.2;Alpha_hung=0.003_0.003_0.01;Tau=0.2_0.8_5;Tau_hung=0.05_0.05_0.15;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    %if N == 8
	%file = '2StrainJuryDelib-3sec;Beta=1;N=8;Binomial_p=0.53;MVM_p=0.6_0.035_0.8;Alpha=0.02_0.02_0.2;Alpha_hung=0.002_0.002_0.01;Tau=0.4_0.6_4;Tau_hung=0.02_0.04_0.15;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    %elseif N == 12
	%file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0.53;MVM_p=0.65_0.04_0.97;Alpha=0.01_0.02_0.2;Alpha_hung=0.002_0.003_0.015;Tau=0.2_0.6_4;Tau_hung=0.02_0.05_0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    %end

    InfluenceCDF = Influencecdf_paramspace(fileDt240,dt,Tboundaries);
    [~, lMVM_240] = bInfluencefit(h, Tboundaries,InfluenceCDF,fileDt240);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%    Comparisons of Distributions   %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [R_MVM_15_60,p_MVM_15_60] = blrtest_l1l2(lMVM_15,lMVM_60,0);
    [R_MVM_15_240,p_MVM_15_240] = blrtest_l1l2(lMVM_15,lMVM_240,0);
    [R_MVM_60_240,p_MVM_60_240] = blrtest_l1l2(lMVM_60,lMVM_240,0);



    R = [R_MVM_15_60 R_MVM_15_240 R_MVM_60_240];
    p = [p_MVM_15_60 p_MVM_15_240 p_MVM_60_240];


end
