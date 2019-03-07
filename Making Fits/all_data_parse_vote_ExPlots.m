function [CDFOR6,CDFOR12,CDFCA6,CDFCA8,CDFCA12] = all_data_parse_vote_ExPlots()
directory = '/export/data/ccbdata/keith/DelibWork/DelibData/';
load(strcat(directory,'NewBin_VoteBinnedData.mat'));
load(strcat(directory,'Influence_Null_5_17_16_1minT_NewCA.mat'));% CHANGE
mindat = 23;

count = 0;

low = 1;
high = 2;

%{
%OR 6
for i=1:1
    count = count + 1;

    data = eval(strcat('VoteBinnedOR6data',num2str(i)));

    b = (0:length(data(1,:)))./60;
    h = data;
    N = numel(h(:,1))-1;
    file = '2StrainJuryDelib-3sec;Beta=1;N=6;Binomial_p=0.5_0.01_0.55;MVM_p=0.86_0.02_1.0;Alpha=0.004_0.003_0.042;Alpha_hung=0_1_0;Tau=2_4_30;HungRatio=0.01_0.07_0.36;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    dt = 60; % dt*N seconds per timestep
    InfluenceCDF = Influencecdf_paramspace(file,dt,b);
    CDFOR6 = GetInfluenceCDF(h,b,InfluenceCDF,file);
    %Params = reshape(all_param_quantiles(count,1:6,2:3),6,2);
    %CDFOR6_Low  =  FindCDF_Influence(h,InfluenceCDF,Params(1,low), Params(2,low),Params(3,low), Params(4,low),Params(5,low), Params(6,low),file);
    %CDFOR6_High =  FindCDF_Influence(h,InfluenceCDF,Params(1,high), Params(2,high),Params(3,high), Params(4,high),Params(5,high), Params(6,high),file);

    save(strcat(directory,'OR6InfluenceVoteTimeBestFitModel_1minT.mat'),'CDFOR6');
end

%OR 12
for i=1:1
    count = count + 1;
    
    data = eval(strcat('VoteBinnedOR12data',num2str(i)));
    
    b = (0:length(data(1,:)))./60;
    h = data;
    N = numel(h(:,1))-1;
    file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0.5_0.01_0.55;MVM_p=0.65_0.04_0.95;Alpha=0.004_0.005_0.03;Alpha_hung=0_1_0;Tau=1_4_20;HungRatio=0.01_0.04_0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
    dt = 60; % dt*N seconds per timestep
    InfluenceCDF = Influencecdf_paramspace(file,dt,b);
    CDFOR12 = GetInfluenceCDF(h,b,InfluenceCDF,file);
    %Params = reshape(all_param_quantiles(count,1:6,2:3),6,2);
    %CDFOR12_Low  =  FindCDF_Influence(h,InfluenceCDF,Params(1,low), Params(2,low),Params(3,low), Params(4,low),Params(5,low), Params(6,low),file);
    %CDFOR12_High =  FindCDF_Influence(h,InfluenceCDF,Params(1,high), Params(2,high),Params(3,high), Params(4,high),Params(5,high), Params(6,high),file);

    save(strcat(directory,'OR12InfluenceVoteTimeBestFitModel_1minT.mat'),'CDFOR12');
end

%CA 6
file = '2StrainJuryDelib-3sec;Beta=1;N=6;Binomial_p=0.5_0.01_0.55;MVM_p=0.86_0.02_1.0;Alpha=0.004_0.003_0.042;Alpha_hung=0_1_0;Tau=2_4_30;HungRatio=0.01_0.07_0.36;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
dt = 60; % dt*N seconds per timestep
Tboundaries = 0:numel(VoteBinnedCA6data(1,1,:));
n = 6;
InfluenceCDF = Influencecdf_paramspace(file,dt,Tboundaries);

    data = zeros(n+1,numel(VoteBinnedCA6data(1,1,:)));
    for i=1:9
        count = count + 1;
    	data = data + reshape(VoteBinnedCA6data(i,:,:),n+1,numel(VoteBinnedCA6data(i,1,:)));
    end    
    if sum(reshape(data,numel(data),1)) > mindat
        b = 0:length(data(1,:));
        No12Hr = b(1:end-1)==12;
        h = data;
        h(:,No12Hr) = 0; 
        CDFCA6 = GetInfluenceCDF(h,b,InfluenceCDF,file);
        %Params = reshape(all_param_quantiles(count,1:6,2:3),6,2);
        %CDFCA6_Low  =  FindCDF_Influence(h,InfluenceCDF,Params(1,low), Params(2,low),Params(3,low), Params(4,low),Params(5,low), Params(6,low),file);
        %CDFCA6_High =  FindCDF_Influence(h,InfluenceCDF,Params(1,high), Params(2,high),Params(3,high), Params(4,high),Params(5,high), Params(6,high),file);
    end
    save(strcat(directory,'CA6InfluenceVoteTimeBestFitModel_1minT.mat'),'CDFCA6');


%CA8
file = '2StrainJuryDelib-3sec;Beta=1;N=8;Binomial_p=0.5_0.01_0.55;MVM_p=0.55_0.03_0.82;Alpha=0.004_0.006_0.05;Alpha_hung=0_1_0;Tau=1_3_25;HungRatio=0.01_0.03_0.2;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
dt = 60; % dt*N seconds per timestep
Tboundaries = 0:numel(VoteBinnedCA8data(1,1,:));
n = 8;
InfluenceCDF = Influencecdf_paramspace(file,dt,Tboundaries);
CDFCA8 = zeros(9,n+1,numel(Tboundaries));
CDFCA8_Low = zeros(9,n+1,numel(Tboundaries));
CDFCA8_High = zeros(9,n+1,numel(Tboundaries));

for i=1:9
    count = count + 1;
    data = reshape(VoteBinnedCA8data(i,:,:),n+1,numel(VoteBinnedCA8data(i,1,:)));
    
    if sum(reshape(data,numel(data),1)) > mindat
        b = (0:length(data(1,:)));
        No12Hr = b(1:end-1)==12;
        h = data;
        h(:,No12Hr) = 0;
        CDFCA8(i,:,:) = GetInfluenceCDF(h,b,InfluenceCDF,file);
        %Params = reshape(all_param_quantiles(count,1:6,2:3),6,2);
        %CDFCA8_Low(i,:,:)  =  FindCDF_Influence(h,InfluenceCDF,Params(1,low), Params(2,low),Params(3,low), Params(4,low),Params(5,low), Params(6,low),file);
        %CDFCA8_High(i,:,:) =  FindCDF_Influence(h,InfluenceCDF,Params(1,high), Params(2,high),Params(3,high), Params(4,high),Params(5,high), Params(6,high),file);
    end
end
save(strcat(directory,'CA8InfluenceVoteTimeBestFitModel_1minT.mat'),'CDFCA8','CDFCA8_Low','CDFCA8_High');
%}
%CA12
file = '2StrainJuryDelib-3sec;Beta=1;N=12;Binomial_p=0.5_0.01_0.55;MVM_p=0.65_0.04_0.95;Alpha=0.004_0.005_0.03;Alpha_hung=0_1_0;Tau=1_4_20;HungRatio=0.01_0.04_0.24;NumRuns=5000;NumTrials=32;CompleteGraph;AllInfected;Trial1.dat';
dt = 60; % dt*N seconds per timestep
Tboundaries = 0:numel(VoteBinnedCA12data(1,1,:));
InfluenceCDF = Influencecdf_paramspace(file,dt,Tboundaries);
n=12;
CDFCA12 = zeros(9,n+1,numel(Tboundaries));
CDFCA12_Low = zeros(9,n+1,numel(Tboundaries));
CDFCA12_High = zeros(9,n+1,numel(Tboundaries));

for i=1:9
    count = count + 1;
    data = reshape(VoteBinnedCA12data(i,:,:),n+1,numel(VoteBinnedCA12data(i,1,:)));
    
    if sum(reshape(data,numel(data),1)) > mindat
        b = 0:length(data(1,:));
        No12Hr = b(1:end-1)==12;
        h = data;
        h(:,No12Hr) = 0;
        CDFCA12(i,:,:) = GetInfluenceCDF(h,b,InfluenceCDF,file);
        %Params = reshape(all_param_quantiles(count,1:6,2:3),6,2);
        %CDFCA12_Low(i,:,:)  =  FindCDF_Influence(h,InfluenceCDF,Params(1,low), Params(2,low),Params(3,low), Params(4,low),Params(5,low), Params(6,low),file);
        %CDFCA12_High(i,:,:) =  FindCDF_Influence(h,InfluenceCDF,Params(1,high), Params(2,high),Params(3,high), Params(4,high),Params(5,high), Params(6,high),file);
    end
end
save(strcat(directory,'CA12InfluenceVoteTimeBestFitModel_1minT.mat'),'CDFCA12','CDFCA12_Low','CDFCA12_High');

save(strcat(directory,'Influence_BestFitModels_1minT.mat'),'CDFOR6','CDFOR12','CDFCA6','CDFCA8','CDFCA12','CDFOR6_Low','CDFOR12_Low','CDFCA6_Low','CDFCA8_Low','CDFCA12_Low','CDFOR6_High','CDFOR12_High','CDFCA6_High','CDFCA8_High','CDFCA12_High');
end
