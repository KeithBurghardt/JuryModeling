(* ::Package:: *)

function [all_p, all_param _est,all_param _errors,all_param _quantiles,all_N,all_Bmin,all_Bmin _quantiles] = all_data _parse _vote _NewInfluence _RatioBinom()
% directory = '/Users/keithburghardt/Google Drive/DelibWork/DelibData/';
directory = '/Users/networklab/Google Drive/DelibWork/DelibData/';

% load(strcat(directory,'NewBin_VoteBinnedData.mat'));
load(strcat(directory,'VoteBinnedData_AllStates _CrimCivil.mat'));
N=12;

all_p = zeros(2+3*9,3+ N+1);
all_param _est = zeros(2+3*9,4,N+1);
all_param _quantiles = zeros((2+3*9),(N+1 + 3 + 6), 4);
all_param _errors = zeros((2+3*9),4,N+1);
all_N = zeros((2+3*9),2,N + 1);
all_Bmin = zeros(2+3*9,N+1);
all_Bmin _quantiles = zeros(2+3*9,N+1,4);
count = 0;
% bootstrap reps
reps = 1000;

mindat = 23;

FullModel = 1;
      NoF = 0;
 NoFAlpha = 0;
    NoFMu = 0;
     NoMu = 0;
  NoQTime = 0;
% else: p = 0.5


StrDescription = '';
if FullModel
    StrDescription = 'FullModel_Q0=0.3';
    file=''
elseif NoF
    StrDescription = 'NoF';
    file=''
elseif NoFAlpha
    StrDescription = 'NoFAlpha';
    file=''
elseif NoFMu
    StrDescription = 'NoFMu';
    file=''
elseif NoMu
    StrDescription = 'NoMu';
    file=''
elseif NoQTime
    StrDescription = 'NoQTime';
    file=''
else
    StrDescription = 'p=0p5';    
    file=''
end



% OR6 civil
for i=1:1
    count = count + 1;

    data = eval('OR6Civil');
    
    b = (0:length(data(1,:)))./60;
    h = data;
    N = numel(h(:,1))-1;
    % file = '2StrainJuryDelib-3sec;Beta=1;N=6;Binomial_p=0.5_ 0.01_ 0.55;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.004_ 0.003_ 0.042;Alpha_hung=0_ 1_ 0;Tau=2_ 4_ 30;HungRatio=0.01_ 0.07_ 0.36;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
	    if FullModel
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_ 1_ 0.5814;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.022;Alpha_hung=0_ 1_ 0;Tau=1_ 2_ 15;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Civil.dat';
        
        % no F
	    elseif NoF
	    	% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6108_ 1_ 0.6108;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.02_ 0.005_ 0.04;Alpha_hung=0_ 1_ 0;Tau=1_ 3_ 16;HungRatio=0.0_ 0.01_ 0.0;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;OR6BinomInitCondit;NoF;Trial1.dat';
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_ 1_ 0.5814;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.022;Alpha_hung=0_ 1_ 0;Tau=1_ 2_ 15;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Civil _NoHungCondit.dat';
	    % no F Alpha
	    elseif NoFAlpha
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_ 1_ 0.5814;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.022;Alpha_hung=0_ 1_ 0;Tau=1_ 2_ 15;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Civil _NoStopHungCondit.dat';
	    % no F Mu
	    elseif NoFMu
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_ 1_ 0.5814;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.022;Alpha_hung=0_ 1_ 0;Tau=1_ 2_ 15;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Civil _NoStubbornHungCondit.dat';
	    elseif NoMu
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_ 1_ 0.5814;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.022;Alpha_hung=0_ 1_ 0;Tau=0_ 2_ 0;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Civil.dat';
        % no q time dependence
        elseif NoQTime
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_ 1_ 0.5814;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.022;Alpha_hung=0_ 1_ 0;Tau=1_ 2_ 15;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;OR6_Civil _NoVoteDependence.dat';
	    % p = 0.5
	    else
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_ 1_ 0.5814;MVM_p=0.5_ 0.02_ 0.5;Alpha=0.010_ 0.002_ 0.022;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Civil.dat';
	    end

    dt = 60; % dt*N seconds per timestep
    InfluenceCDF = Influencecdf_paramspace(file,dt,b);
    [Params, BestL] = bInfluencefit(h,b,InfluenceCDF,file);
    disp('BEST FIT:');
    disp(Params);
    disp('BEST L 10^-4:');
    test=BestL;
    test(BestL==-100000000) = log(10^-4);
    disp(sum(sum(test)));

    disp('BEST L 10^-7:');
    test=BestL;
    test(BestL==-100000000) = log(10^-7);
    disp(sum(sum(test)));

    disp('BEST L 10^-11:');
    test=BestL;
    test(BestL==-100000000) = log(10^-11);
    disp(sum(sum(test)));

    disp('BEST L 10^-14:');
    test=BestL;
    test(BestL==-100000000) = log(10^-14);
    disp(sum(sum(test)));

    
    [p, param_est,param_errors,param_quantiles,N_tail,Bmin,Bmin_quantiles] = bfinddist_vote _NEW(h, b,dt,file,InfluenceCDF,reps);
    save(strcat(directory,'OR6VoteTimeDistModel_Civil _ 1minT_',StrDescription,'_New _ 1000reps.mat'),'p','param_est','param_errors','param_quantiles','N_tail','Bmin','Bmin_quantiles');
    all_p(count,1:N+1+3) = p;
    all_param _est(count,:,1:N+1) = param_est;
    all_param _errors(count,:,1:+N+1) = param_errors;
    all_param _quantiles(count,1:(N+1 + 3 + 6),:) = param_quantiles;
    all_N(count,:,1:N+1) = N_tail;
    all_Bmin(count,1:N+1) = Bmin;
    all_Bmin _quantiles(count,1:N+1,:) = Bmin_quantiles;
    
end


% OR6 criminal
for i=1:1
    count = count + 1;

    data = eval('OR6Criminal');
    
    b = (0:length(data(1,:)))./60;
    h = data;
    N = numel(h(:,1))-1;
    % file = '2StrainJuryDelib-3sec;Beta=1;N=6;Binomial_p=0.5_ 0.01_ 0.55;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.004_ 0.003_ 0.042;Alpha_hung=0_ 1_ 0;Tau=2_ 4_ 30;HungRatio=0.01_ 0.07_ 0.36;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
	    if FullModel
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.015_ 0.002_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.002_ 0.01;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;OR6_Criminal.dat';
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.025_ 0.002_ 0.035;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.002_ 0.01;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;OR6_Criminal.dat';
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.025_ 0.002_ 0.035;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;OR6_Criminal.dat';
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.02_ 0.002_ 0.03;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;OR6_Criminal.dat';
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.02_ 0.002_ 0.03;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal.dat';
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.015_ 0.002_ 0.027;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal.dat';
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.024;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal.dat';
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.024;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal.dat';
         % no F
	    elseif NoF
            file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.024;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal _NoHungCondit.dat';

        % no F Alpha
	    elseif NoFAlpha
            file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.024;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal _NoStopHungCondit.dat';
	    % no F Mu
	    elseif NoFMu
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.024;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal _NoStubbornHungCondit.dat';    
        % Mu = 0
        elseif NoMu
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal.dat';
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.010_ 0.002_ 0.024;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;OR6_Criminal _NoVoteDependence.dat';
	    % p = 0.5
        else
            file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_ 1_ 0.6390;MVM_p=0.5_ 0.02_ 0.5;Alpha=0.010_ 0.002_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal.dat';
        end

    dt = 60; % dt*N seconds per timestep
    InfluenceCDF = Influencecdf_paramspace(file,dt,b);
    [Params, BestL] = bInfluencefit(h,b,InfluenceCDF,file);
    disp('BEST FIT:');
    disp(Params);
    disp('BEST L 10^-4:');
    test=BestL;
    test(BestL==-100000000) = log(10^-4);
    disp(sum(sum(test)));

    disp('BEST L 10^-7:');
    test=BestL;
    test(BestL==-100000000) = log(10^-7);
    disp(sum(sum(test)));

    disp('BEST L 10^-11:');
    test=BestL;
    test(BestL==-100000000) = log(10^-11);
    disp(sum(sum(test)));

    disp('BEST L 10^-14:');
    test=BestL;
    test(BestL==-100000000) = log(10^-14);
    disp(sum(sum(test)));

    
    [p, param_est,param_errors,param_quantiles,N_tail,Bmin,Bmin_quantiles] = bfinddist_vote _NEW(h, b,dt,file,InfluenceCDF,reps);
    save(strcat(directory,'OR6VoteTimeDistModel_Criminal _ 1minT_',StrDescription,'_New _ 1000reps.mat'),'p','param_est','param_errors','param_quantiles','N_tail','Bmin','Bmin_quantiles');
    all_p(count,1:N+1+3) = p;
    all_param _est(count,:,1:N+1) = param_est;
    all_param _errors(count,:,1:+N+1) = param_errors;
    all_param _quantiles(count,1:(N+1 + 3 + 6),:) = param_quantiles;
    all_N(count,:,1:N+1) = N_tail;
    all_Bmin(count,1:N+1) = Bmin;
    all_Bmin _quantiles(count,1:N+1,:) = Bmin_quantiles;
    
end


% OR 12 civil
for i=1:1
    count = count + 1;
    
    data = eval('OR12Civil');
    
    b = (0:length(data(1,:)))./60;
    h = data;
    N = numel(h(:,1))-1;
    % file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0.5_ 0.01_ 0.55;MVM_p=0.65_ 0.04_ 0.95;Alpha=0.004_ 0.005_ 0.03;Alpha_hung=0_ 1_ 0;Tau=1_ 4_ 20;HungRatio=0.01_ 0.04_ 0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
	    if FullModel
    	    file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.89_ 0.005_ 0.93;Alpha=0.015_ 0.002_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.002_ 0.01;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;OR12_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.015_ 0.002_ 0.025;Alpha_hung=0_ 1_ 0;Tau=5_ 10_ 55;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;OR12_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.015_ 0.002_ 0.025;Alpha_hung=0_ 1_ 0;Tau=5_ 10_ 55;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 7_ 40;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil.dat';
	    % no F
	    elseif NoF
        	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5452_ 1_ 0.5452;MVM_p=0.85_ 0.02_ 0.95;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.001_ 0.02_ 0.001;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12BinomInitCondit;NoF.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil _NoHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=50_ 10_ 101;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil _NoHungCondit.dat';
	    % no F Alpha
	    elseif NoFAlpha
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5452_ 1_ 0.5452;MVM_p=0.85_ 0.02_ 0.95;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.001_ 0.02_ 0.101;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12BinomInitCondit;NoFAlpha.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil _NoStopHungCondit.dat';
	    % no F Mu
	    elseif NoFMu
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5452_ 1_ 0.5452;MVM_p=0.85_ 0.02_ 0.95;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.001_ 0.02_ 0.101;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12BinomInitCondit;NoFMu.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil _NoStubbornHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=50_ 10_ 101;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil _NoStubbornHungCondit.dat';
	    elseif NoMu
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5452_ 1_ 0.5452;MVM_p=0.85_ 0.02_ 0.97;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.001_ 0.02_ 0.101;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.88_ 0.02_ 0.98;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0_ 10_ 0;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil.dat';
        
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5452_ 1_ 0.5452;MVM_p=0.89_ 0.005_ 0.92;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=1;CompleteGraph;AllInfected;OR12_BinomInitCondit _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;OR12_Civil _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.84_ 0.01_ 0.91;Alpha=0.002_ 0.002_ 0.010;Alpha_hung=0_ 1_ 0;Tau=40_ 10_ 101;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;OR12_Civil _NoVoteDependence.dat';
	    % p = 0.5
	    else
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5452_ 1_ 0.5452;MVM_p=0.5_ 0.02_ 0.5;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.001_ 0.02_ 0.101;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_ 1_ 0.5169;MVM_p=0.5_ 0.01_ 0.5;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil.dat';
	    end

    dt = 60; % dt*N seconds per timestep
    InfluenceCDF = Influencecdf_paramspace(file,dt,b);
    [Params, BestL] = bInfluencefit(h,b,InfluenceCDF,file);
    disp('BEST FIT:');
    disp(Params);
    disp('BEST L 10^-4:');
    test=BestL;
    test(BestL==-100000000) = log(10^-4);
    disp(sum(sum(test)));

    disp('BEST L 10^-7:');
    test=BestL;
    test(BestL==-100000000) = log(10^-7);
    disp(sum(sum(test)));

    disp('BEST L 10^-11:');
    test=BestL;
    test(BestL==-100000000) = log(10^-11);
    disp(sum(sum(test)));

    disp('BEST L 10^-14:');
    test=BestL;
    test(BestL==-100000000) = log(10^-14);
    disp(sum(sum(test)));

    
    [p, param_est,param_errors,param_quantiles,N_tail,Bmin,Bmin_quantiles] = bfinddist_vote _NEW(h, b,dt,file,InfluenceCDF,reps);
    save(strcat(directory,'OR12VoteTimeDistModel_Civil _ 1minT_',StrDescription,'_New _ 1000reps.mat'),'p','param_est','param_errors','param_quantiles','N_tail','Bmin','Bmin_quantiles');
    all_p(count,1:N+1+3) = p;
    all_param _est(count,:,1:N+1) = param_est;
    all_param _errors(count,:,1:+N+1) = param_errors;
    all_param _quantiles(count,1:(N+1 + 3 + 6),:) = param_quantiles;
    all_N(count,:,1:N+1) = N_tail;
    all_Bmin(count,1:N+1) = Bmin;
    all_Bmin _quantiles(count,1:N+1,:) = Bmin_quantiles;
    
end


% OR 12 criminal
for i=1:1
    count = count + 1;
    
    data = eval('OR12Criminal');
    
    b = (0:length(data(1,:)))./60;
    h = data;
    N = numel(h(:,1))-1;
    % file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0.5_ 0.01_ 0.55;MVM_p=0.65_ 0.04_ 0.95;Alpha=0.004_ 0.005_ 0.03;Alpha_hung=0_ 1_ 0;Tau=1_ 4_ 20;HungRatio=0.01_ 0.04_ 0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
	    if FullModel
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_ 1_ 0.5904;MVM_p=0.85_ 0.02_ 0.95;Alpha=0.004_ 0.002_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 20_ 151;HungRatio=0.0_ 0.2_ 1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal.dat';
	    % no F
	    elseif NoF
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_ 1_ 0.5904;MVM_p=0.85_ 0.02_ 0.95;Alpha=0.004_ 0.002_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 2_ 15;HungRatio=0.0_ 0.05_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal _NoHungCondit.dat';
	    % no F Alpha
	    elseif NoFAlpha
               file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_ 1_ 0.5904;MVM_p=0.85_ 0.02_ 0.95;Alpha=0.004_ 0.002_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 2_ 15;HungRatio=0.0_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal _NoStopHungCondit.dat';
	    % no F Mu
	    elseif NoFMu
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_ 1_ 0.5904;MVM_p=0.85_ 0.02_ 0.95;Alpha=0.004_ 0.002_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 2_ 15;HungRatio=0.0_ 0.05_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal _NoStubbornHungCondit.dat';
            
        % Mu = 0
	    elseif NoMu
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_ 1_ 0.5904;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.004_ 0.002_ 0.016;Alpha_hung=0_ 1_ 0;Tau=0_ 2_ 0;HungRatio=0.0_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_ 1_ 0.5904;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.004_ 0.002_ 0.016;Alpha_hung=0_ 1_ 0;Tau=0_ 2_ 0;HungRatio=0.0_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal.dat';
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_ 1_ 0.5904;MVM_p=0.85_ 0.02_ 0.95;Alpha=0.004_ 0.002_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 5_ 40;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;OR12_Criminal _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_ 1_ 0.5904;MVM_p=0.85_ 0.02_ 0.95;Alpha=0.004_ 0.002_ 0.016;Alpha_hung=0_ 1_ 0;Tau=20_ 20_ 100;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;OR12_Criminal _NoVoteDependence.dat';
	    % p = 0.5
        else
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_ 1_ 0.5904;MVM_p=0.5_ 0.02_ 0.5;Alpha=0.004_ 0.002_ 0.016;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal.dat';
	    end

    dt = 60; % dt*N seconds per timestep
    InfluenceCDF = Influencecdf_paramspace(file,dt,b);
    [Params, BestL] = bInfluencefit(h,b,InfluenceCDF,file);
    disp('BEST FIT:');
    disp(Params);
    disp('BEST L 10^-4:');
    test=BestL;
    test(BestL==-100000000) = log(10^-4);
    disp(sum(sum(test)));

    disp('BEST L 10^-7:');
    test=BestL;
    test(BestL==-100000000) = log(10^-7);
    disp(sum(sum(test)));

    disp('BEST L 10^-11:');
    test=BestL;
    test(BestL==-100000000) = log(10^-11);
    disp(sum(sum(test)));

    disp('BEST L 10^-14:');
    test=BestL;
    test(BestL==-100000000) = log(10^-14);
    disp(sum(sum(test)));

    
    [p, param_est,param_errors,param_quantiles,N_tail,Bmin,Bmin_quantiles] = bfinddist_vote _NEW(h, b,dt,file,InfluenceCDF,reps);
    save(strcat(directory,'OR12VoteTimeDistModel_Criminal _ 1minT_',StrDescription,'_New4 _ 1000reps.mat'),'p','param_est','param_errors','param_quantiles','N_tail','Bmin','Bmin_quantiles');
    all_p(count,1:N+1+3) = p;
    all_param _est(count,:,1:N+1) = param_est;
    all_param _errors(count,:,1:+N+1) = param_errors;
    all_param _quantiles(count,1:(N+1 + 3 + 6),:) = param_quantiles;
    all_N(count,:,1:N+1) = N_tail;
    all_Bmin(count,1:N+1) = Bmin;
    all_Bmin _quantiles(count,1:N+1,:) = Bmin_quantiles;
    
end

disp(all_p(1:count,:));
disp(all_param _est(1:count,:,:));
disp(all_param _errors(1:count,:,:));
disp(all_param _quantiles(1:count,:,:));
disp(all_N(1:count,:,:));
disp(all_Bmin(1:count,:));
disp(all_Bmin _quantiles(1:count,:,:));

% CA 6
	    % file = '2StrainJuryDelib-3sec;Beta=1;N=6;Binomial_p=0.5_ 0.01_ 0.55;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.004_ 0.003_ 0.042;Alpha_hung=0_ 1_ 0;Tau=2_ 4_ 30;HungRatio=0.01_ 0.07_ 0.36;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
            % Full Model
            if FullModel
		% file='2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0_ 1_ 0;MVM_p=0.86_ 0.02_ 1.0;Alpha=0.004_ 0.003_ 0.042;Alpha_hung=0_ 1_ 0;Tau=2_ 4_ 30;HungRatio=0.01_ 0.07_ 0.36;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.005_ 0.002_ 0.015;Alpha_hung=0_ 1_ 0;Tau=4_ 10_ 44;HungRatio=0.001_ 0.1_ 0.81;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.002_ 0.001_ 0.008;Alpha_hung=0_ 1_ 0;Tau=4_ 10_ 44;HungRatio=0.5_ 0.1_ 1.0;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.002_ 0.001_ 0.008;Alpha_hung=0_ 1_ 0;Tau=4_ 10_ 44;HungRatio=0.2_ 0.1_ 0.8;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.008_ 0.003_ 0.032;Alpha_hung=0_ 1_ 0;Tau=4_ 10_ 44;HungRatio=0.2_ 0.1_ 1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.002_ 0.003_ 0.014;Alpha_hung=0_ 1_ 0;Tau=4_ 10_ 44;HungRatio=0.2_ 0.1_ 1.0;NumRuns=5000;NumTrials=32;Q0=0.5;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.008_ 0.003_ 0.032;Alpha_hung=0_ 1_ 0;Tau=4_ 10_ 44;HungRatio=0.2_ 0.1_ 1.0;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.002_ 0.003_ 0.014;Alpha_hung=0_ 1_ 0;Tau=4_ 10_ 44;HungRatio=0.2_ 0.1_ 1.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=4_ 10_ 44;HungRatio=0.2_ 0.1_ 1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
            % No f
            elseif NoF
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.90_ 0.025_ 1.00;Alpha=0.004_ 0.003_ 0.016;Alpha_hung=0_ 1_ 0;Tau=0.5_ 7_ 40;HungRatio=0.0_ 0.04_ 0.0;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoF;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.92_ 0.02_ 1.00;Alpha=0.001_ 0.003_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.0_ 0.04_ 0.0;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoF;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=10_ 10_ 70;HungRatio=1_ 0.1_ 1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;NoF.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=10_ 10_ 70;HungRatio=0_ 0.1_ 0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;NoF.dat';
            % No f for Alpha
            elseif NoFAlpha
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.90_ 0.025_ 1.00;Alpha=0.004_ 0.003_ 0.016;Alpha_hung=0_ 1_ 0;Tau=0.5_ 7_ 40;HungRatio=0.01_ 0.09_ 0.48;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFAlpha;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.85_ 0.025_ 1.00;Alpha=0.001_ 0.003_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.1_ 0.1_ 0.7;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFAlpha;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.85_ 0.025_ 1.00;Alpha=0.001_ 0.003_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.1_ 0.1_ 0.7;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFAlpha;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.85_ 0.025_ 1.00;Alpha=0.001_ 0.003_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.1_ 0.1_ 0.7;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFAlpha;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.75_ 0.025_ 0.95;Alpha=0.008_ 0.004_ 0.024;Alpha_hung=0_ 1_ 0;Tau=1_ 15_ 100;HungRatio=0.0_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFAlpha;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.75_ 0.025_ 0.95;Alpha=0.008_ 0.004_ 0.024;Alpha_hung=0_ 1_ 0;Tau=1_ 15_ 100;HungRatio=0.0_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFAlpha;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.004_ 0.004_ 0.024;Alpha_hung=0_ 1_ 0;Tau=50_ 15_ 150;HungRatio=0.0_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFAlpha;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.004_ 0.004_ 0.024;Alpha_hung=0_ 1_ 0;Tau=50_ 15_ 150;HungRatio=0.0_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFAlpha;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=10_ 10_ 70;HungRatio=0.2_ 0.2_ 1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFAlpha.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.2_ 0.2_ 1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFAlpha.dat';
            % No f for stubbornness
            elseif NoFMu
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.90_ 0.025_ 1.00;Alpha=0.004_ 0.003_ 0.016;Alpha_hung=0_ 1_ 0;Tau=0.5_ 7_ 40;HungRatio=0.01_ 0.09_ 0.48;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.85_ 0.025_ 1.00;Alpha=0.001_ 0.003_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.1_ 0.1_ 0.7;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.85_ 0.025_ 1.00;Alpha=0.001_ 0.003_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.1_ 0.1_ 0.7;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.85_ 0.025_ 1.00;Alpha=0.001_ 0.003_ 0.016;Alpha_hung=0_ 1_ 0;Tau=1_ 10_ 61;HungRatio=0.1_ 0.1_ 0.7;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.75_ 0.025_ 0.95;Alpha=0.008_ 0.004_ 0.024;Alpha_hung=0_ 1_ 0;Tau=30_ 10_ 100;HungRatio=0.0_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.75_ 0.025_ 0.95;Alpha=0.008_ 0.004_ 0.024;Alpha_hung=0_ 1_ 0;Tau=30_ 10_ 100;HungRatio=0.0_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.75_ 0.025_ 0.95;Alpha=0.008_ 0.004_ 0.024;Alpha_hung=0_ 1_ 0;Tau=30_ 10_ 100;HungRatio=0.0_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.003_ 0.003_ 0.018;Alpha_hung=0_ 1_ 0;Tau=30_ 10_ 100;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.003_ 0.003_ 0.018;Alpha_hung=0_ 1_ 0;Tau=30_ 10_ 100;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.003_ 0.003_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 5_ 31;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=10_ 10_ 70;HungRatio=0.2_ 0.2_ 1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.2_ 0.2_ 1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;NoFMu.dat';
	    elseif NoMu
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.01_ 1.0;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.2_ 0.2_ 1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit.dat';
        % no q time dependence
        elseif NoQTime
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.007_ 1.0;Alpha=0.007_ 0.002_ 0.016;Alpha_hung=0_ 1_ 0;Tau=10_ 5_ 40;HungRatio=0.5_ 0.1_ 1.0;NumRuns=5000;NumTrials=32;Q0=1;CompleteGraph;AllInfected;CA6_BinomInitCondit.dat';
        file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.95_ 0.007_ 1.0;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=10_ 10_ 80;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1;CompleteGraph;AllInfected;CA6_BinomInitCondit _NoVoteDependence.dat';
            % p = 0.5
            else
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.5_ 0.03_ 0.5;Alpha=0.03_ 0.01_ 0.09;Alpha_hung=0_ 1_ 0;Tau=0.1_ 0.2_ 1;HungRatio=0.001_ 0.004_ 0.005;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.5_ 0.025_ 0.5;Alpha=0.01_ 0.01_ 0.06;Alpha_hung=0_ 1_ 0;Tau=1_ 5_ 31;HungRatio=0.01_ 0.05_ 0.31;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		% file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.5_ 0.025_ 0.5;Alpha=0.04_ 0.01_ 0.1;Alpha_hung=0_ 1_ 0;Tau=10_ 10_ 61;HungRatio=0.0_ 0.01_ 0.05;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.5_ 0.025_ 0.5;Alpha=0.04_ 0.01_ 0.1;Alpha_hung=0_ 1_ 0;Tau=10_ 10_ 61;HungRatio=0.0_ 0.01_ 0.05;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.5_ 0.01_ 0.5;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.2_ 0.2_ 1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit.dat';
		file = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_ 1_ 0.4545;MVM_p=0.5_ 0.01_ 0.5;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.01_ 0.05_ 0.41;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit.dat';
            end

        dt = 60; % dt*N seconds per timestep
        Tboundaries = 0:numel(VoteBinnedCA6data(1,1,:));
        InfluenceCDF = Influencecdf_paramspace(file,dt,Tboundaries);

    n = 6;
    data = zeros(n+1,numel(VoteBinnedCA6data(1,1,:)));
    for i=1:9
        count = count + 1;
    	data = data + reshape(VoteBinnedCA6data(i,:,:),n+1,numel(VoteBinnedCA6data(i,1,:)));
    end    
    if sum(reshape(data,numel(data),1)) > mindat
        b = 0:length(data(1,:));
        
        h = data;
        % remove hung juries
        % h(floor(n/4)+2:n-floor(n/4),:) = 0;
        N = numel(h(:,1))-1;
	[Params, BestL] = bInfluencefit(h,Tboundaries,InfluenceCDF,file);
	disp('BEST FIT:');
	disp(Params);

    	disp('BEST L 10^-4:');
    	test=BestL;
    	test(BestL==-100000000) = log(10^-4);
    	disp(sum(sum(test)));

    	disp('BEST L 10^-7:');
    	test=BestL;
    	test(BestL==-100000000) = log(10^-7);
    	disp(sum(sum(test)));

    	disp('BEST L 10^-11:');
    	test=BestL;
    	test(BestL==-100000000) = log(10^-11);
    	disp(sum(sum(test)));

    	disp('BEST L 10^-14:');
    	test=BestL;
        test(BestL==-100000000) = log(10^-14);
    	disp(sum(sum(test)));


   
          	
        [p, param_est,param_errors,param_quantiles,N_tail,Bmin,Bmin_quantiles] = bfinddist_vote _NEW(h, b,dt,file,InfluenceCDF,reps);
	all_p(count,1:N+1+3) = p;
    	all_param _est(count,:,1:N+1) = param_est;
    	all_param _errors(count,:,1:+N+1) = param_errors;
    	all_param _quantiles(count,1:(N+1 + 3 + 6),:) = param_quantiles;
    	all_N(count,:,1:N+1) = N_tail;
    	all_Bmin(count,1:N+1) = Bmin;
    	all_Bmin _quantiles(count,1:N+1,:) = Bmin_quantiles;
	
    end

save(strcat(directory,'CA6VoteTimeDistModel_FBinom _ 1minT_',StrDescription,'_New.mat'),'p','param_est','param_errors','param_quantiles','all_N','Bmin','Bmin_quantiles');

disp(all_p(1:count,:));
disp(all_param _est(1:count,:,:));
disp(all_param _errors(1:count,:,:));
disp(all_param _quantiles(1:count,:,:));
disp(all_N(1:count,:,:));
disp(all_Bmin(1:count,:));
disp(all_Bmin _quantiles(1:count,:,:));



% CA8
% file = '2StrainJuryDelib-3sec;Beta=1;N=8;Binomial_p=0.5_ 0.01_ 0.55;MVM_p=0.55_ 0.03_ 0.82;Alpha=0.004_ 0.006_ 0.05;Alpha_hung=0_ 1_ 0;Tau=1_ 3_ 25;HungRatio=0.01_ 0.03_ 0.2;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';

dt = 60; % dt*N seconds per timestep
Tboundaries = 0:numel(VoteBinnedCA8data(1,1,:));

for i=4:4
    count = count + 1;
    n = 8;
    data = reshape(VoteBinnedCA8data(i,:,:),n+1,numel(VoteBinnedCA8data(i,1,:)));
    if sum(reshape(data,numel(data),1)) > mindat
        b = (0:length(data(1,:)));
        h = data;
        N = numel(h(:,1))-1;
	if i == 4
   	    % file = '2StrainJuryDelib-3sec;Beta=1;N=8;Binomial_p=0_ 1_ 0;MVM_p=0.55_ 0.03_ 0.82;Alpha=0.004_ 0.006_ 0.05;Alpha_hung=0_ 1_ 0;Tau=1_ 6_ 40;HungRatio=0.01_ 0.03_ 0.2;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA8_ 1BinomInitCondit;Trial1.dat';
	    % full model
	    if FullModel
   	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5583_ 1_ 0.5583;MVM_p=0.75_ 0.04_ 0.99;Alpha=0.01_ 0.006_ 0.05;Alpha_hung=0_ 1_ 0;Tau=0.2_ 0.2_ 2;HungRatio=0.001_ 0.05_ 0.41;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_BinomInitCondit;Trial1.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.78_ 0.02_ 0.88;Alpha=0.01_ 0.002_ 0.02;Alpha_hung=0_ 1_ 0;Tau=1_ 0.1_ 2;HungRatio=0.0_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.78_ 0.02_ 0.88;Alpha=0.015_ 0.005_ 0.035;Alpha_hung=0_ 1_ 0;Tau=1_ 0.1_ 2;HungRatio=0.1_ 0.1_ 1.0;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.78_ 0.02_ 0.88;Alpha=0.015_ 0.005_ 0.035;Alpha_hung=0_ 1_ 0;Tau=1_ 0.1_ 2;HungRatio=0.0_ 0.1_ 0.4;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.78_ 0.02_ 0.88;Alpha=0.015_ 0.005_ 0.035;Alpha_hung=0_ 1_ 0;Tau=1_ 0.1_ 2;HungRatio=0.0_ 0.1_ 0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.70_ 0.05_ 0.95;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.5_ 0.2_ 2;HungRatio=0.0_ 0.1_ 0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.70_ 0.05_ 0.95;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.0_ 0.1_ 0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.87_ 0.03_ 1.0;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.0_ 0.1_ 0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.82_ 0.03_ 0.97;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.0_ 0.1_ 0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.75_ 0.05_ 1.0;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.0_ 0.1_ 0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
        % no F
	    elseif NoF
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5583_ 1_ 0.5583;MVM_p=0.88_ 0.02_ 1.0;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.001_ 0.04_ 0.001;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_BinomInitCondit;NoF.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.70_ 0.05_ 0.95;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil _NoHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=4_ 2_ 16;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil _NoHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.82_ 0.03_ 1.0;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=1_ 2_ 11;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil _NoHungCondit.dat';
	    % no F Alpha
	    elseif NoFAlpha
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5583_ 1_ 0.5583;MVM_p=0.8_ 0.04_ 1.0;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.001_ 0.04_ 0.201;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_BinomInitCondit;NoFAlpha.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.70_ 0.05_ 0.95;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.0_ 0.1_ 0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil _NoStopHungCondit.dat';
	    % no F Mu
	    elseif NoFMu
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5583_ 1_ 0.5583;MVM_p=0.8_ 0.04_ 1.0;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.001_ 0.04_ 0.201;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_BinomInitCondit;NoFMu.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.70_ 0.05_ 0.95;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil _NoStubbornHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.90_ 0.02_ 1.0;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=4_ 2_ 15;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil _NoStubbornHungCondit.dat';
	    elseif NoMu
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5583_ 1_ 0.5583;MVM_p=0.8_ 0.04_ 1.0;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=0_ 0.3_ 0;HungRatio=0.001_ 0.04_ 0.201;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.75_ 0.05_ 1.0;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.0_ 0.1_ 0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5583_ 1_ 0.5583;MVM_p=0.88_ 0.02_ 1.0;Alpha=0.005_ 0.002_ 0.015;Alpha_hung=0_ 1_ 0;Tau=0.5_ 0.1_ 1;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=1;CA8_ 1_BinomInitCondit _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.87_ 0.03_ 1.0;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA8_ 1_Civil _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.002_ 0.002_ 0.012;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA8_ 1_Civil _NoVoteDependence.dat';
	    % p = 0.5
	    else
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5583_ 1_ 0.5583;MVM_p=0.5_ 0.04_ 0.5;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.001_ 0.04_ 0.201;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_ 1_ 0.5667;MVM_p=0.5_ 0.03_ 0.5;Alpha=0.008_ 0.003_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.1_ 0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 1_Civil.dat';
	    end

	elseif i == 5
   	    % file = '2StrainJuryDelib-3sec;Beta=1;N=8;Binomial_p=0_ 1_ 0;MVM_p=0.55_ 0.03_ 0.82;Alpha=0.004_ 0.006_ 0.05;Alpha_hung=0_ 1_ 0;Tau=1_ 6_ 40;HungRatio=0.01_ 0.03_ 0.2;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA8_ 2BinomInitCondit;Trial1.dat';
	    % full model
	    if FullModel
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5960_ 1_ 0.5960;MVM_p=0.69_ 0.04_ 0.93;Alpha=0.005_ 0.005_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 0.5_ 6;HungRatio=0.04_ 0.05_ 0.31;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_BinomInitCondit;Trial1.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.78_ 0.02_ 0.88;Alpha=0.005_ 0.002_ 0.015;Alpha_hung=0_ 1_ 0;Tau=1_ 0.1_ 2;HungRatio=0.1_ 0.02_ 0.2;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA8_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.78_ 0.02_ 0.88;Alpha=0.013_ 0.002_ 0.023;Alpha_hung=0_ 1_ 0;Tau=1_ 0.1_ 2;HungRatio=0.1_ 0.02_ 0.2;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA8_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.78_ 0.02_ 0.88;Alpha=0.013_ 0.002_ 0.023;Alpha_hung=0_ 1_ 0;Tau=1_ 0.1_ 2;HungRatio=0.1_ 0.02_ 0.2;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.75_ 0.05_ 0.95;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 8;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.75_ 0.05_ 1.0;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 6;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.75_ 0.05_ 1.0;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 6;HungRatio=0.0_ 0.2_ 1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_Civil.dat';
	    % no F
	    elseif NoF
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5960_ 1_ 0.5960;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.005_ 0.003_ 0.02;Alpha_hung=0_ 1_ 0;Tau=0.5_ 0.5_ 3;HungRatio=0.05_ 0.05_ 0.05;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_BinomInitCondit;NoF.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.75_ 0.05_ 1.0;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_Civil _NoHungCondit.dat';
	    % no F Alpha
	    elseif NoFAlpha
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5960_ 1_ 0.5960;MVM_p=0.7_ 0.03_ 0.9;Alpha=0.005_ 0.003_ 0.02;Alpha_hung=0_ 1_ 0;Tau=0.5_ 0.5_ 3;HungRatio=0.05_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_BinomInitCondit;NoFAlpha.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.75_ 0.05_ 1.0;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 6;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_Civil _NoStopHungCondit.dat';
	    % no F Mu
	    elseif NoFMu
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5960_ 1_ 0.5960;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.005_ 0.003_ 0.02;Alpha_hung=0_ 1_ 0;Tau=0.5_ 0.5_ 3;HungRatio=0.05_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_BinomInitCondit;NoFMu.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.75_ 0.05_ 1.0;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_Civil _NoStubbornHungCondit.dat';

	    elseif NoMu
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5960_ 1_ 0.5960;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.005_ 0.003_ 0.02;Alpha_hung=0_ 1_ 0;Tau=0.0_ 0.5_ 0;HungRatio=0.05_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.75_ 0.05_ 1.0;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0_ 0.7_ 0;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_Civil.dat';
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5960_ 1_ 0.5960;MVM_p=0.88_ 0.02_ 1.0;Alpha=0.005_ 0.002_ 0.015;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=1;CA8_ 2_BinomInitCondit _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.75_ 0.05_ 1.0;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 6;HungRatio=0.0_ 0.2_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA8_ 2_Civil _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.75_ 0.05_ 1.0;Alpha=0.002_ 0.002_ 0.012;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 6;HungRatio=0.0_ 0.2_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA8_ 2_Civil _NoVoteDependence.dat';
	    % p = 0.5
	    else
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5960_ 1_ 0.5960;MVM_p=0.5_ 0.03_ 0.5;Alpha=0.005_ 0.003_ 0.02;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.001_ 0.03_ 0.201;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_ 1_ 0.6003;MVM_p=0.5_ 0.05_ 0.5;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.2_ 1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 2_Civil.dat';
	    end

	elseif i == 6
   	    % file = '2StrainJuryDelib-3sec;Beta=1;N=8;Binomial_p=0_ 1_ 0;MVM_p=0.55_ 0.03_ 0.82;Alpha=0.004_ 0.006_ 0.05;Alpha_hung=0_ 1_ 0;Tau=1_ 6_ 40;HungRatio=0.01_ 0.03_ 0.2;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA8_ 3BinomInitCondit;Trial1.dat';
	    % full model
	    if FullModel
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_ 1_ 0.5579;MVM_p=0.60_ 0.05_ 0.9;Alpha=0.002_ 0.003_ 0.015;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 9;HungRatio=0.001_ 0.02_ 0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 3_BinomInitCondit;Trial1.dat';
	    % no F
	    elseif NoF
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_ 1_ 0.5579;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.001_ 0.02_ 0.001;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 3_BinomInitCondit;NoF.dat';
	    % no F Alpha
	    elseif NoFAlpha
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_ 1_ 0.5579;MVM_p=0.7_ 0.03_ 0.9;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0.5_ 0.4_ 4;HungRatio=0.001_ 0.02_ 0.101;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 3_BinomInitCondit;NoFAlpha.dat';
	    % no F Mu
	    elseif NoFMu
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_ 1_ 0.5579;MVM_p=0.82_ 0.03_ 1.0;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0.5_ 0.4_ 4;HungRatio=0.001_ 0.02_ 0.101;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 3_BinomInitCondit;NoFMu.dat';
	    elseif NoMu
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_ 1_ 0.5579;MVM_p=0.85_ 0.03_ 1;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0_ 0.4_ 0;HungRatio=0.001_ 0.02_ 0.101;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 3_BinomInitCondit.dat';
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_ 1_ 0.5579;MVM_p=0.88_ 0.02_ 1.0;Alpha=0.002_ 0.002_ 0.012;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 2;HungRatio=0.0_ 0.01_ 0.0;NumRuns=5000;NumTrials=32;Q0=1;CA8_ 3_BinomInitCondit _NoVoteDependence.dat';
	    % p = 0.5
	    else
	    	file = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_ 1_ 0.5579;MVM_p=0.5_ 0.03_ 0.5;Alpha=0.001_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=0.1_ 0.1_ 1;HungRatio=0.001_ 0.02_ 0.101;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_ 3_BinomInitCondit.dat';
           
	    end
    end
        InfluenceCDF = Influencecdf_paramspace(file,dt,Tboundaries);
		
    	[Params, BestL] = bInfluencefit(h,Tboundaries,InfluenceCDF,file);
    	disp('BEST FIT:');
        disp(Params);

    	disp('BEST L 10^-4:');
    	test=BestL;
    	test(BestL==-100000000) = log(10^-4);
    	disp(sum(sum(test)));

    	disp('BEST L 10^-7:');
    	test=BestL;
    	test(BestL==-100000000) = log(10^-7);
    	disp(sum(sum(test)));

    	disp('BEST L 10^-11:');
    	test=BestL;
    	test(BestL==-100000000) = log(10^-11);
    	disp(sum(sum(test)));

    	disp('BEST L 10^-14:');
    	test=BestL;
        test(BestL==-100000000) = log(10^-14);
    	disp(sum(sum(test)));

        
        
        [p, param_est,param_errors,param_quantiles,N_tail,Bmin,Bmin_quantiles] = bfinddist_vote _NEW(h, b,dt,file,InfluenceCDF,reps);
        all_p(count,1:N+1+3) = p;
    	all_param _est(count,:,1:N+1) = param_est;
    	all_param _errors(count,:,1:+N+1) = param_errors;
    	all_param _quantiles(count,1:(N+1 + 3 + 6),:) = param_quantiles;
    	all_N(count,:,1:N+1) = N_tail;
    	all_Bmin(count,1:N+1) = Bmin;
    	all_Bmin _quantiles(count,1:N+1,:) = Bmin_quantiles;
        
    end
end

save(strcat(directory,'CA8Vote_',StrDescription,'_New4 _ 1000reps-4.mat'),'all_p','all_param _est','all_param _errors','all_param _quantiles','all_N','all_Bmin','all_Bmin _quantiles');

disp(all_p(1:count,:));
disp(all_param _est(1:count,:,:));
disp(all_param _errors(1:count,:,:));
disp(all_param _quantiles(1:count,:,:));
disp(all_N(1:count,:,:));
disp(all_Bmin(1:count,:));
disp(all_Bmin _quantiles(1:count,:,:));



% CA12
dt = 60; % dt*N seconds per timestep
Tboundaries = 0:numel(VoteBinnedCA12data(1,1,:));

for i = 4:4
    % if i~=5
    count = count + 1;
    n = 12;
    data = reshape(VoteBinnedCA12data(i,:,:),n+1,numel(VoteBinnedCA12data(i,1,:)));
    
    if sum(reshape(data,numel(data),1)) > mindat
        b = 0:length(data(1,:));
        h = data;
        N = numel(h(:,1))-1;
	if i == 4
   	    % file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0_ 1_ 0;MVM_p=0.65_ 0.04_ 0.95;Alpha=0.004_ 0.005_ 0.03;Alpha_hung=0_ 1_ 0;Tau=1_ 4_ 20;HungRatio=0.01_ 0.04_ 0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA12_ 1BinomInitCondit;Trial1.dat';
	    % Full Model
	    if FullModel
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.75_ 0.05_ 0.95;Alpha=0.012_ 0.002_ 0.03;Alpha_hung=0_ 1_ 0;Tau=0.2_ 0.2_ 2;HungRatio=0.001_ 0.05_ 0.31;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_BinomInitCondit;Trial1.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.80_ 0.02_ 0.90;Alpha=0.012_ 0.001_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 0.1_ 2;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.80_ 0.02_ 0.90;Alpha=0.015_ 0.002_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 0.1_ 2;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.80_ 0.02_ 0.90;Alpha=0.015_ 0.002_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 0.1_ 2;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.80_ 0.02_ 0.90;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 5;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.75_ 0.03_ 0.90;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.79_ 0.04_ 0.99;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil.dat';
            % no F
	    elseif NoF
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.88_ 0.02_ 1.0;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_BinomInitCondit;NoF.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.75_ 0.03_ 0.90;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil _NoHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil _NoHungCondit.dat';
	    elseif NoFAlpha
	    % No F Alpha
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.88_ 0.02_ 1.0;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.1_ 0.1_ 1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_BinomInitCondit;NoFAlpha.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.75_ 0.03_ 0.90;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil _NoStopHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.75_ 0.03_ 0.90;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil _NoStopHungCondit.dat';
	    % No F Mu
	    elseif NoFMu
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.88_ 0.02_ 1.0;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.001_ 0.1_ 0.501;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_BinomInitCondit;NoFMu.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.75_ 0.03_ 0.90;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil _NoStubbornHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil _NoStubbornHungCondit.dat';
	    elseif NoMu
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.88_ 0.02_ 1.0;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=0_ 0.3_ 0;HungRatio=0.1_ 0.1_ 1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.80_ 0.03_ 0.95;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0_ 0.7_ 0;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0_ 0.7_ 0;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil.dat';
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.88_ 0.02_ 1.0;Alpha=0.002_ 0.002_ 0.012;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.03_ 0.0;NumRuns=5000;NumTrials=32;Q0=1;CA12_ 1_BinomInitCondit _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.79_ 0.04_ 0.99;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_ 1_Civil _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.79_ 0.04_ 0.99;Alpha=0.002_ 0.002_ 0.012;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_ 1_Civil _NoVoteDependence.dat';
	    % p = 0.5
	    else
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.5_ 0.02_ 0.5;Alpha=0.005_ 0.004_ 0.025;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.01_ 0.05_ 0.41;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.5_ 0.04_ 0.5;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_ 1_ 0.5379;MVM_p=0.5_ 0.04_ 0.5;Alpha=0.014_ 0.002_ 0.024;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 1_Civil.dat';
            
	    end
	elseif i == 5
   	    % file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0_ 1_ 0;MVM_p=0.65_ 0.04_ 0.95;Alpha=0.004_ 0.005_ 0.03;Alpha_hung=0_ 1_ 0;Tau=1_ 4_ 20;HungRatio=0.01_ 0.04_ 0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA12_ 2BinomInitCondit;Trial1.dat';
	    % Full Model
	    if FullModel
 	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5192_ 1_ 0.5192;MVM_p=0.75_ 0.02_ 0.91;Alpha=0.005_ 0.002_ 0.015;Alpha_hung=0_ 1_ 0;Tau=0.5_ 0.5_ 7;HungRatio=0.05_ 0.05_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_BinomInitCondit;Trial1.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.80_ 0.02_ 0.90;Alpha=0.009_ 0.001_ 0.013;Alpha_hung=0_ 1_ 0;Tau=1_ 0.3_ 3;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.80_ 0.02_ 0.90;Alpha=0.01_ 0.002_ 0.02;Alpha_hung=0_ 1_ 0;Tau=1_ 0.3_ 3;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.74_ 0.02_ 0.84;Alpha=0.01_ 0.002_ 0.02;Alpha_hung=0_ 1_ 0;Tau=1_ 0.3_ 3;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.74_ 0.02_ 0.84;Alpha=0.01_ 0.002_ 0.02;Alpha_hung=0_ 1_ 0;Tau=1_ 0.3_ 3;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.74_ 0.04_ 0.90;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 6;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.74_ 0.04_ 0.90;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 6;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil.dat';
	    % No F
	    elseif NoF
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5192_ 1_ 0.5192;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.005_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.01_ 0.05_ 0.01;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_BinomInitCondit;NoF.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.74_ 0.04_ 0.90;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil _NoHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil _NoHungCondit.dat';
	    % No F Alpha
	    elseif NoFAlpha
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5192_ 1_ 0.5192;MVM_p=0.7_ 0.03_ 0.9;Alpha=0.005_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.01_ 0.1_ 0.71;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_BinomInitCondit;NoFAlpha.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.74_ 0.04_ 0.90;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 6;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil _NoStopHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.74_ 0.04_ 0.90;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 6;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil _NoStopHungCondit.dat';
	    % No F Mu
	    elseif NoFMu
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5192_ 1_ 0.5192;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.005_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.01_ 0.05_ 0.31;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_BinomInitCondit;NoFMu.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.74_ 0.04_ 0.90;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil _NoStubbornHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil _NoStubbornHungCondit.dat';
	    elseif NoMu
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5192_ 1_ 0.5192;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.005_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.01_ 0.1_ 0.71;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.78_ 0.04_ 0.94;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil.dat';
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5192_ 1_ 0.5192;MVM_p=0.88_ 0.02_ 1.0;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.05_ 0.0;NumRuns=5000;NumTrials=32;Q0=1;CA12_ 2_BinomInitCondit _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.74_ 0.04_ 0.90;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 6;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_ 2_Civil _NoVoteDependence.dat';
	    % p = 0.5
	    else
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5192_ 1_ 0.5192;MVM_p=0.5_ 0.03_ 0.5;Alpha=0.005_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.01_ 0.1_ 0.71;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_ 1_ 0.5329;MVM_p=0.5_ 0.04_ 0.5;Alpha=0.005_ 0.002_ 0.013;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 2_Civil.dat';
	    end
	elseif i == 6
   	    % file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0_ 1_ 0;MVM_p=0.65_ 0.04_ 0.95;Alpha=0.004_ 0.005_ 0.03;Alpha_hung=0_ 1_ 0;Tau=1_ 4_ 20;HungRatio=0.01_ 0.04_ 0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA12_ 3BinomInitCondit;Trial1.dat';
	    if FullModel
	    	% file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5381_ 1_ 0.5381;MVM_p=0.70_ 0.03_ 0.93;Alpha=0.005_ 0.003_ 0.02;Alpha_hung=0_ 1_ 0;Tau=0.2_ 0.5_ 5;HungRatio=0.001_ 0.05_ 0.41;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_BinomInitCondit;Trial1.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.75_ 0.02_ 0.85;Alpha=0.005_ 0.001_ 0.011;Alpha_hung=0_ 1_ 0;Tau=1_ 0.3_ 3;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 3_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.75_ 0.02_ 0.85;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 0.3_ 3;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 3_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.75_ 0.02_ 0.85;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=1_ 0.3_ 3;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.7_ 0.05_ 0.9;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 8;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.7_ 0.05_ 0.9;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 8;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.7_ 0.05_ 0.9;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil.dat';
	    % No F
	    elseif NoF
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5381_ 1_ 0.5381;MVM_p=0.74_ 0.04_ 0.95;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=1_ 4_ 21;HungRatio=0.05_ 0.05_ 0.05;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_BinomInitCondit;NoF.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.7_ 0.05_ 0.9;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 8;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil _NoHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 8;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil _NoHungCondit.dat';
	    % No F Alpha
	    elseif NoFAlpha
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5381_ 1_ 0.5381;MVM_p=0.74_ 0.04_ 0.95;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 7;HungRatio=0.01_ 0.05_ 0.31;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_BinomInitCondit;NoFAlpha.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.7_ 0.05_ 0.9;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 8;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil _NoStopHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.7_ 0.05_ 0.9;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 8;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil _NoStopHungCondit.dat';
	    % No F Mu
	    elseif NoFMu
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5381_ 1_ 0.5381;MVM_p=0.74_ 0.04_ 0.95;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=1_ 4_ 21;HungRatio=0.01_ 0.05_ 0.31;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_BinomInitCondit;NoFMu.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.7_ 0.05_ 0.9;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 8;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil _NoStubbornHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 8;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil _NoStubbornHungCondit.dat';
	    elseif NoMu
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5381_ 1_ 0.5381;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.01_ 0.05_ 0.31;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.75_ 0.05_ 0.95;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil.dat';
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5381_ 1_ 0.5381;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=1_ 0.3_ 3;HungRatio=0.0_ 0.03_ 0.0;NumRuns=5000;NumTrials=32;Q0=1;CA12_ 3_BinomInitCondit _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.7_ 0.05_ 0.9;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_ 3_Civil _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.002_ 0.002_ 0.010;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.7_ 4;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_ 3_Civil _NoVoteDependence.dat';
	    % p = 0.5
	    else
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5381_ 1_ 0.5381;MVM_p=0.5_ 0.04_ 0.5;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.01_ 0.05_ 0.31;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_ 1_ 0.5437;MVM_p=0.5_ 0.05_ 0.5;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.1_ 0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 3_Civil.dat';
	    end

	elseif i == 7
   	    % file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0_ 1_ 0;MVM_p=0.65_ 0.04_ 0.95;Alpha=0.004_ 0.005_ 0.03;Alpha_hung=0_ 1_ 0;Tau=1_ 4_ 20;HungRatio=0.01_ 0.04_ 0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA12_ 4BinomInitCondit;Trial1.dat';
	    % Full Model
	    if FullModel
	    	% file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.63_ 0.03_ 0.81;Alpha=0.005_ 0.003_ 0.026;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 6;HungRatio=0.001_ 0.03_ 0.21;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_BinomInitCondit;Trial1.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.02_ 0.78;Alpha=0.005_ 0.001_ 0.011;Alpha_hung=0_ 1_ 0;Tau=2_ 1_ 8;HungRatio=0.0_ 0.02_ 0.12;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 4_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.02_ 0.78;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=2_ 1_ 8;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 4_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.02_ 0.78;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=5_ 3_ 20;HungRatio=0.1_ 0.05_ 0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 4_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.02_ 0.78;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=5_ 3_ 20;HungRatio=0.0_ 0.05_ 0.2;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_ 4_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.02_ 0.78;Alpha=0.008_ 0.002_ 0.018;Alpha_hung=0_ 1_ 0;Tau=5_ 3_ 20;HungRatio=0.0_ 0.05_ 0.2;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.03_ 0.86;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=1_ 2_ 12;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.03_ 0.86;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 9;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_Civil.dat';
	    % No F
	    elseif NoF
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.75_ 0.04_ 0.95;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.05_ 0.05_ 0.05;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_BinomInitCondit;NoF.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.03_ 0.86;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 9;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_Civil _NoHungCondit.dat';
	    % No F Alpha
	    elseif NoFAlpha
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.65_ 0.04_ 0.85;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.01_ 0.04_ 0.21;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_BinomInitCondit;NoFAlpha.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.03_ 0.86;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 9;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_Civil _NoStopHungCondit.dat';
	    % No F Mu
	    elseif NoFMu
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.75_ 0.04_ 0.95;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.01_ 0.04_ 0.21;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_BinomInitCondit;NoFMu.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.03_ 0.86;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 9;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_Civil _NoStubbornHungCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 9;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_Civil _NoStubbornHungCondit.dat';
	    elseif NoMu
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.01_ 0.04_ 0.21;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.75_ 0.03_ 0.93;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_Civil.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.88_ 0.03_ 1.0;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_Civil.dat';
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.8_ 0.02_ 0.92;Alpha=0.001_ 0.001_ 0.006;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.0_ 0.02_ 0.0;NumRuns=5000;NumTrials=32;Q0=1;CA12_ 4_BinomInitCondit _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.68_ 0.03_ 0.86;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 9;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_ 4_Civil _NoVoteDependence.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.82_ 0.03_ 1.0;Alpha=0.001_ 0.001_ 0.007;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 9;HungRatio=0.0_ 0.1_ 0.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_ 4_Civil _NoVoteDependence.dat';
	    % p = 0.5
	    else
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5318_ 1_ 0.5318;MVM_p=0.5_ 0.04_ 0.5;Alpha=0.003_ 0.003_ 0.021;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.01_ 0.04_ 0.21;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_BinomInitCondit.dat';
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_ 1_ 0.5456;MVM_p=0.5_ 0.03_ 0.5;Alpha=0.003_ 0.002_ 0.011;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.0_ 0.1_ 0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 4_Civil.dat';
	    end

	elseif i == 8
   	    % file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0_ 1_ 0;MVM_p=0.65_ 0.04_ 0.95;Alpha=0.004_ 0.005_ 0.03;Alpha_hung=0_ 1_ 0;Tau=1_ 4_ 20;HungRatio=0.01_ 0.04_ 0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;CA12_ 5BinomInitCondit;Trial1.dat';
	    % Full Model
	    if FullModel
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_ 1_ 0.5265;MVM_p=0.63_ 0.03_ 0.84;Alpha=0.002_ 0.001_ 0.007;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 8;HungRatio=0.001_ 0.03_ 0.21;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 5_BinomInitCondit;Trial1.dat';

	    % No f
	    elseif NoF
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_ 1_ 0.5265;MVM_p=0.83_ 0.03_ 1.0;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.001_ 0.05_ 0.001;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 5_BinomInitCondit;NoF.dat';
	    % No f for Alpha
	    elseif NoFAlpha
    	    file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_ 1_ 0.5265;MVM_p=0.7_ 0.03_ 0.9;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.001_ 0.05_ 0.301;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 5_BinomInitCondit;NoFAlpha.dat';
	    % No f for stubbornness
	    elseif NoFMu
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_ 1_ 0.5265;MVM_p=0.83_ 0.03_ 1.0;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=1_ 1_ 10;HungRatio=0.001_ 0.05_ 0.301;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 5_BinomInitCondit;NoFMu.dat';
	    elseif NoMu
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_ 1_ 0.5265;MVM_p=0.85_ 0.03_ 1.0;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=0_ 1_ 0;HungRatio=0.001_ 0.05_ 0.301;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 5_BinomInitCondit.dat';
        % no q time dependence
        elseif NoQTime
            file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_ 1_ 0.5265;MVM_p=0.9_ 0.02_ 1.0;Alpha=0.0005_ 0.0005_ 0.004;Alpha_hung=0_ 1_ 0;Tau=0.5_ 1_ 6;HungRatio=0.0_ 0.005_ 0.0;NumRuns=5000;NumTrials=32;Q0=1;CA12_ 5_BinomInitCondit _NoVoteDependence.dat';
	    % p = 0.5
	    else
	    	file = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_ 1_ 0.5265;MVM_p=0.5_ 0.03_ 0.5;Alpha=0.002_ 0.002_ 0.01;Alpha_hung=0_ 1_ 0;Tau=0.3_ 0.3_ 3;HungRatio=0.001_ 0.05_ 0.301;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_ 5_BinomInitCondit.dat';
	    end
	end
	InfluenceCDF = Influencecdf_paramspace(file,dt,Tboundaries);
    	[Params, BestL] = bInfluencefit(h,Tboundaries,InfluenceCDF,file);
    	disp('BEST FIT:');
    	disp(Params);

    	disp('BEST L 10^-4:');
    	test=BestL;
    	test(BestL==-100000000) = log(10^-4);
    	disp(sum(sum(test)));

    	disp('BEST L 10^-7:');
    	test=BestL;
    	test(BestL==-100000000) = log(10^-7);
    	disp(sum(sum(test)));

    	disp('BEST L 10^-11:');
    	test=BestL;
    	test(BestL==-100000000) = log(10^-11);
    	disp(sum(sum(test)));

    	disp('BEST L 10^-14:');
    	test=BestL;
        test(BestL==-100000000) = log(10^-14);
    	disp(sum(sum(test)));
		    
        [p, param_est,param_errors,param_quantiles,N_tail,Bmin,Bmin_quantiles] = bfinddist_vote _NEW(h, b,dt,file,InfluenceCDF,reps);
        all_p(count,1:N+1+3) = p;
    	all_param _est(count,:,1:N+1) = param_est;
    	all_param _errors(count,:,1:+N+1) = param_errors;
    	all_param _quantiles(count,1:(N+1 + 3 + 6),:) = param_quantiles;
    	all_N(count,:,1:N+1) = N_tail;
    	all_Bmin(count,1:N+1) = Bmin;
    	all_Bmin _quantiles(count,1:N+1,:) = Bmin_quantiles;
        
    end
    % end
end

% save(strcat(directory,'CA12Vote_','New2_',StrDescription,'_ 1000reps_ 4,6.mat'),'all_p','all_param _est','all_param _errors','all_param _quantiles','all_N','all_Bmin','all_Bmin _quantiles');

disp(all_p(1:count,:));
disp(all_param _est(1:count,:,:));
disp(all_param _errors(1:count,:,:));
disp(all_param _quantiles(1:count,:,:));
disp(all_N(1:count,:,:));
disp(all_Bmin(1:count,:));
disp(all_Bmin _quantiles(1:count,:,:));

% save(strcat(directory,'Influence_FBinom _',StrDescription,'_New.mat'),'all_p','all_param _est','all_param _errors','all_param _quantiles','all_N','all_Bmin','all_Bmin _quantiles');




end



