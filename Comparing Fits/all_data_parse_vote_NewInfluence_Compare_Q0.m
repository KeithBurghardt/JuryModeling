function [all_Rp,all_R] = all_data_parse_vote_NewInfluence_Compare_Q0()
directory = '~/Google Drive/DelibWork/DelibData/';
%load(strcat(directory,'NewBin_VoteBinnedData.mat'));
load(strcat(directory,'VoteBinnedData_AllStates_CrimCivil.mat'));
%load(strcat(directory,'ECHRMatFiles.mat'));

all_Rp = zeros(2 + 9 * 3,3);
all_R = zeros(2 + 9 * 3,3);
count = 0;

mindat = 23;


%OR 6 Civil
for i=1:1
    count = count + 1;

    %data = eval(strcat('VoteBinnedOR6data',num2str(i)));
    data = eval(strcat('OR6Civil'));
    b = (0:length(data(1,:)))./60;
    h = data;
    fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_1_0.5814;MVM_p=0.9_0.02_1.0;Alpha=0.010_0.002_0.022;Alpha_hung=0_1_0;Tau=1_2_15;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;OR6_Civil.dat';
    fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_1_0.5814;MVM_p=0.9_0.02_1.0;Alpha=0.010_0.002_0.022;Alpha_hung=0_1_0;Tau=1_2_15;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Civil.dat';
    fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_1_0.5814;MVM_p=0.9_0.02_1.0;Alpha=0.010_0.002_0.022;Alpha_hung=0_1_0;Tau=1_2_15;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;OR6_Civil.dat';

    %Dts: CHANGE TIMESTEPS in bCompareQ0s
    %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_1_0.5814;MVM_p=0.9_0.02_1.0;Alpha=0.0025_0.0005_0.0055;Alpha_hung=0_1_0;Tau=1_2_15;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Civil.dat';
    %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.5814_1_0.5814;MVM_p=0.9_0.02_1.0;Alpha=0.040_0.008_0.088;Alpha_hung=0_1_0;Tau=1_2_15;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Civil.dat';
    

    [R,Rp] = bCompareQ0s(h, b,fileQ01,fileQ03,fileQ05);
    save(strcat(directory,'OR6Vote_Civil_Q0=0.1,0.3,1.0.mat'),'R','Rp');
    all_R(count,:) = R;
    all_Rp(count,:) = Rp;
end

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));
%}
%{
%OR 6 Criminal
for i=1:1
    count = count + 1;

    %data = eval(strcat('VoteBinnedOR6data',num2str(i)));
    data = eval(strcat('OR6Criminal'));
    b = (0:length(data(1,:)))./60;
    h = data;
    %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_1_0.6390;MVM_p=0.86_0.02_1.0;Alpha=0.010_0.002_0.024;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;OR6_Criminal.dat';
    fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_1_0.6390;MVM_p=0.86_0.02_1.0;Alpha=0.010_0.002_0.024;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal.dat';
    %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_1_0.6390;MVM_p=0.86_0.02_1.0;Alpha=0.010_0.002_0.024;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;OR6_Criminal.dat';

    %Dts: CHANGE TIMESTEPS
    
    fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_1_0.6390;MVM_p=0.86_0.02_1.0;Alpha=0.0025_0.0005_0.0060;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal.dat';
    fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.6390_1_0.6390;MVM_p=0.86_0.02_1.0;Alpha=0.040_0.008_0.096;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR6_Criminal.dat';

    [R,Rp] = bCompareQ0s(h, b,fileQ01,fileQ03,fileQ05);
    save(strcat(directory,'OR6Vote_Criminal_Dt=15,60,240.mat'),'R','Rp');
    all_R(count,:) = R;
    all_Rp(count,:) = Rp;
end

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));
%}
%{
%OR 12 Civil
for i=1:1
    count = count + 1;
    
    data = eval('OR12Civil');
    
    b = (0:length(data(1,:)))./60;
    h = data;
    
    %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_1_0.5169;MVM_p=0.84_0.01_0.91;Alpha=0.015_0.002_0.025;Alpha_hung=0_1_0;Tau=5_10_55;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;OR12_Civil.dat';
    fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_1_0.5169;MVM_p=0.84_0.01_0.91;Alpha=0.008_0.002_0.018;Alpha_hung=0_1_0;Tau=1_10_61;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil.dat';
    %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_1_0.5169;MVM_p=0.84_0.01_0.91;Alpha=0.008_0.002_0.018;Alpha_hung=0_1_0;Tau=1_10_61;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;OR12_Civil.dat';
    %Dts
    
    fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_1_0.5169;MVM_p=0.84_0.01_0.91;Alpha=0.0020_0.0005_0.0045;Alpha_hung=0_1_0;Tau=1_10_61;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil.dat';
    fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5169_1_0.5169;MVM_p=0.84_0.01_0.91;Alpha=0.032_0.008_0.072;Alpha_hung=0_1_0;Tau=1_10_61;HungRatio=0.0_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Civil.dat';
    

    [R,Rp] = bCompareQ0s(h, b,fileQ01,fileQ03,fileQ05);
    save(strcat(directory,'OR12Vote_Civil_Dt=15,60,240.mat'),'R','Rp');
    all_R(count,:) = R;
    all_Rp(count,:) = Rp;
end

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));
%}
%{
%OR 12 Criminal
for i=1:1
    count = count + 1;
    
    data = eval('OR12Criminal');
    
    b = (0:length(data(1,:)))./60;
    h = data;
    %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_1_0.5904;MVM_p=0.85_0.02_0.95;Alpha=0.004_0.002_0.016;Alpha_hung=0_1_0;Tau=1_10_81;HungRatio=0.0_0.2_1.0;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;OR12_Criminal.dat';
    fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_1_0.5904;MVM_p=0.85_0.02_0.95;Alpha=0.004_0.002_0.016;Alpha_hung=0_1_0;Tau=1_20_151;HungRatio=0.0_0.2_1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal.dat';
    %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_1_0.5904;MVM_p=0.85_0.02_0.95;Alpha=0.004_0.002_0.016;Alpha_hung=0_1_0;Tau=1_20_151;HungRatio=0.0_0.2_1.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;OR12_Criminal.dat';
    
    %Dts
    fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_1_0.5904;MVM_p=0.85_0.02_0.95;Alpha=0.001_0.0005_0.004;Alpha_hung=0_1_0;Tau=1_2_15;HungRatio=0.0_0.05_0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal.dat';
    fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_1_0.5904;MVM_p=0.85_0.02_0.95;Alpha=0.016_0.008_0.064;Alpha_hung=0_1_0;Tau=1_2_15;HungRatio=0.0_0.05_0.3;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal.dat';
    
    

    [R,Rp] = bCompareQ0s(h, b,fileQ01,fileQ03,fileQ05);
    save(strcat(directory,'OR12Vote_Criminal_Dt=15,60,240.mat'),'R','Rp');
    all_R(count,:) = R;
    all_Rp(count,:) = Rp;
end

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));
%}
%}
%{
%CA 6
%for i=1:9
%    count = count + 1;
    n = 6;
    data = zeros(n+1,numel(VoteBinnedCA6data(1,1,:)));
    for i=1:9
        count = count + 1;
    	data = data + reshape(VoteBinnedCA6data(i,:,:),n+1,numel(VoteBinnedCA6data(i,1,:)));
    end    
    b = 0:length(data(1,:));    
    h = data;

    fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_1_0.4545;MVM_p=0.95_0.01_1.0;Alpha=0.002_0.003_0.014;Alpha_hung=0_1_0;Tau=4_10_44;HungRatio=0.2_0.1_1.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
    fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_1_0.4545;MVM_p=0.95_0.01_1.0;Alpha=0.008_0.003_0.032;Alpha_hung=0_1_0;Tau=4_10_44;HungRatio=0.2_0.1_1.0;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
    fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_1_0.4545;MVM_p=0.95_0.01_1.0;Alpha=0.002_0.002_0.01;Alpha_hung=0_1_0;Tau=4_10_44;HungRatio=0.2_0.1_1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
    %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_1_0.4545;MVM_p=0.95_0.01_1.0;Alpha=0.002_0.003_0.014;Alpha_hung=0_1_0;Tau=4_10_44;HungRatio=0.2_0.1_1.0;NumRuns=5000;NumTrials=32;Q0=0.5;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';

    %Dts
    %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_1_0.4545;MVM_p=0.95_0.01_1.0;Alpha=0.008_0.008_0.04;Alpha_hung=0_1_0;Tau=4_10_44;HungRatio=0.2_0.1_1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';
    %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=6;Binomial_p=0.4545_1_0.4545;MVM_p=0.95_0.01_1.0;Alpha=0.0005_0.0005_0.0025;Alpha_hung=0_1_0;Tau=4_10_44;HungRatio=0.2_0.1_1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA6BinomInitCondit;Trial1.dat';

    [R,Rp] = bCompareQ0s(h, b,fileQ01,fileQ03,fileQ05);
    all_R(count,:) = R;
    all_Rp(count,:) = Rp;
%end
save(strcat(directory,'CA6VoteTimeDistModel_Q0=0.1,0.3,1.mat'),'all_R','all_Rp');

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));

%}
%}
%CA8
%{
for i=4:5
    count = count + 1;
    n = 8;
    data = reshape(VoteBinnedCA8data(i,:,:),n+1,numel(VoteBinnedCA8data(i,1,:)));
    if i ==4
        %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_1_0.5667;MVM_p=0.78_0.02_0.88;Alpha=0.015_0.005_0.035;Alpha_hung=0_1_0;Tau=1_0.1_2;HungRatio=0.0_0.1_0.4;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA8_1_Civil.dat';
        fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_1_0.5667;MVM_p=0.75_0.05_1.0;Alpha=0.008_0.003_0.024;Alpha_hung=0_1_0;Tau=0.5_1_6;HungRatio=0.0_0.1_0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_1_Civil.dat';
        %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_1_0.5667;MVM_p=0.70_0.05_0.95;Alpha=0.008_0.003_0.024;Alpha_hung=0_1_0;Tau=0.5_1_6;HungRatio=0.0_0.1_0.4;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA8_1_Civil.dat';
        
   
        %Dts
        
    	fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_1_0.5667;MVM_p=0.70_0.05_0.95;Alpha=0.002_0.00075_0.006;Alpha_hung=0_1_0;Tau=0.5_1_6;HungRatio=0.0_0.1_0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_1_Civil.dat';
        fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5667_1_0.5667;MVM_p=0.70_0.05_0.95;Alpha=0.032_0.012_0.096;Alpha_hung=0_1_0;Tau=0.5_1_6;HungRatio=0.0_0.1_0.4;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_1_Civil.dat';
   
    elseif i==5
        %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_1_0.6003;MVM_p=0.78_0.02_0.88;Alpha=0.013_0.002_0.023;Alpha_hung=0_1_0;Tau=1_0.1_2;HungRatio=0.1_0.02_0.2;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA8_2_Civil.dat';
        fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_1_0.6003;MVM_p=0.75_0.05_1.0;Alpha=0.008_0.002_0.018;Alpha_hung=0_1_0;Tau=0.3_0.7_6;HungRatio=0.0_0.2_1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_2_Civil.dat';
        %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_1_0.6003;MVM_p=0.75_0.05_1.0;Alpha=0.008_0.002_0.018;Alpha_hung=0_1_0;Tau=0.3_0.7_6;HungRatio=0.0_0.2_1.0;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA8_2_Civil.dat';
  
        %Dts
        fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_1_0.6003;MVM_p=0.75_0.05_1.0;Alpha=0.002_0.0005_0.0045;Alpha_hung=0_1_0;Tau=0.3_0.7_6;HungRatio=0.0_0.2_1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_2_Civil.dat';
        fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.6003_1_0.6003;MVM_p=0.75_0.05_1.0;Alpha=0.032_0.008_0.072;Alpha_hung=0_1_0;Tau=0.3_0.7_6;HungRatio=0.0_0.2_1.0;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_2_Civil.dat';
    elseif i==6
        fileQ05 = '';%'2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_1_0.5579;MVM_p=0.60_0.05_0.9;Alpha=0.001_0.0015_0.01;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.001_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA8_3_BinomInitCondit;Trial1.dat';
        fileQ01 = '';%'2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_1_0.5579;MVM_p=0.60_0.05_0.9;Alpha=0.005_0.005_0.025;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.001_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA8_3_BinomInitCondit;Trial1.dat';
        fileQ03 = '';%'2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_1_0.5579;MVM_p=0.60_0.05_0.9;Alpha=0.002_0.003_0.015;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.001_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_3_BinomInitCondit;Trial1.dat';
        %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_1_0.5579;MVM_p=0.60_0.05_0.9;Alpha=0.001_0.0015_0.01;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.001_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.5;CompleteGraph;AllInfected;CA8_3_BinomInitCondit;Trial1.dat';

        %Dts
        %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_1_0.5579;MVM_p=0.60_0.05_0.9;Alpha=0.008_0.012_0.06;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.001_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_3_BinomInitCondit;Trial1.dat';
        %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=8;Binomial_p=0.5579_1_0.5579;MVM_p=0.60_0.05_0.9;Alpha=0.0005_0.00075_0.0037;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.001_0.02_0.1;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA8_3_BinomInitCondit;Trial1.dat';

    end
    if sum(reshape(data,numel(data),1)) > mindat
        b = (0:length(data(1,:)));
        h = data;
        [R,Rp] = bCompareQ0s(h, b,fileQ01,fileQ03,fileQ05);
    	all_R(count,:) = R;
    	all_Rp(count,:) = Rp;
    end
end

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));
save(strcat(directory,'CA8VoteTime_New_Dt=15,60,240_4-5.mat'),'all_R','all_Rp');
%}
%{
%CA12
for i=5:5
    
    count = count + 1;
    n = 12;
    data = reshape(VoteBinnedCA12data(i,:,:),n+1,numel(VoteBinnedCA12data(i,1,:)));
    if i ==4
        %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_1_0.5379;MVM_p=0.80_0.02_0.90;Alpha=0.015_0.002_0.025;Alpha_hung=0_1_0;Tau=1_0.1_2;HungRatio=0.1_0.05_0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_1_Civil.dat';
        fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_1_0.5379;MVM_p=0.79_0.04_0.99;Alpha=0.008_0.002_0.018;Alpha_hung=0_1_0;Tau=0.3_0.7_4;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_1_Civil.dat';
        %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_1_0.5379;MVM_p=0.75_0.03_0.90;Alpha=0.008_0.002_0.018;Alpha_hung=0_1_0;Tau=0.3_0.7_4;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_1_Civil.dat';

        %Dts
        
        fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_1_0.5379;MVM_p=0.79_0.04_0.99;Alpha=0.002_0.0005_0.0045;Alpha_hung=0_1_0;Tau=0.3_0.7_4;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_1_Civil.dat';
        fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5379_1_0.5379;MVM_p=0.79_0.04_0.99;Alpha=0.032_0.008_0.072;Alpha_hung=0_1_0;Tau=0.3_0.7_4;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_1_Civil.dat';
        
    elseif i==5
        %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_1_0.5329;MVM_p=0.74_0.02_0.84;Alpha=0.01_0.002_0.02;Alpha_hung=0_1_0;Tau=1_0.3_3;HungRatio=0.1_0.05_0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_2_Civil.dat';
        fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_1_0.5329;MVM_p=0.74_0.04_0.90;Alpha=0.005_0.002_0.013;Alpha_hung=0_1_0;Tau=1_1_6;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_2_Civil.dat';
        %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_1_0.5329;MVM_p=0.74_0.04_0.90;Alpha=0.005_0.002_0.013;Alpha_hung=0_1_0;Tau=1_1_6;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_2_Civil.dat';

        %Dts
        fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_1_0.5329;MVM_p=0.74_0.04_0.90;Alpha=0.00125_0.00050_0.00325;Alpha_hung=0_1_0;Tau=1_1_6;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_2_Civil.dat';
        fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5329_1_0.5329;MVM_p=0.74_0.04_0.90;Alpha=0.020_0.008_0.052;Alpha_hung=0_1_0;Tau=1_1_6;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_2_Civil.dat';
        
    elseif i==6
        %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_1_0.5437;MVM_p=0.75_0.02_0.85;Alpha=0.008_0.002_0.018;Alpha_hung=0_1_0;Tau=1_0.3_3;HungRatio=0.1_0.05_0.3;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_3_Civil.dat';
        fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_1_0.5437;MVM_p=0.7_0.05_0.9;Alpha=0.003_0.002_0.011;Alpha_hung=0_1_0;Tau=0.3_0.7_4;HungRatio=0.0_0.1_0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_3_Civil.dat';
        %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_1_0.5437;MVM_p=0.7_0.05_0.9;Alpha=0.003_0.002_0.011;Alpha_hung=0_1_0;Tau=1_1_8;HungRatio=0.0_0.1_0.5;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_3_Civil.dat';

        %Dts
        
        fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_1_0.5437;MVM_p=0.7_0.05_0.9;Alpha=0.00075_0.00050_0.00275;Alpha_hung=0_1_0;Tau=1_1_8;HungRatio=0.0_0.1_0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_3_Civil.dat';
        fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5437_1_0.5437;MVM_p=0.7_0.05_0.9;Alpha=0.012_0.008_0.044;Alpha_hung=0_1_0;Tau=1_1_8;HungRatio=0.0_0.1_0.5;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_3_Civil.dat';
        
    elseif i==7
        %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_1_0.5456;MVM_p=0.68_0.02_0.78;Alpha=0.008_0.002_0.018;Alpha_hung=0_1_0;Tau=5_3_20;HungRatio=0.0_0.05_0.2;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_4_Civil.dat';
        fileQ03 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_1_0.5456;MVM_p=0.68_0.03_0.86;Alpha=0.003_0.002_0.011;Alpha_hung=0_1_0;Tau=0.5_1_9;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_4_Civil.dat';
        %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_1_0.5456;MVM_p=0.68_0.03_0.86;Alpha=0.003_0.002_0.011;Alpha_hung=0_1_0;Tau=0.5_1_9;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_4_Civil.dat';

        %Dts
        
        fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_1_0.5456;MVM_p=0.68_0.03_0.86;Alpha=0.00075_0.00050_0.00275;Alpha_hung=0_1_0;Tau=0.5_1_9;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_4_Civil.dat';
        fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5456_1_0.5456;MVM_p=0.68_0.03_0.86;Alpha=0.012_0.008_0.044;Alpha_hung=0_1_0;Tau=0.5_1_9;HungRatio=0.0_0.1_0.6;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_4_Civil.dat';
        
    elseif i==8
        fileQ05 = '';%'2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_1_0.5265;MVM_p=0.63_0.03_0.84;Alpha=0.001_0.001_0.007;Alpha_hung=0_1_0;Tau=0.5_0.5_5;HungRatio=0.001_0.03_0.21;NumRuns=5000;NumTrials=32;Q0=1.0;CompleteGraph;AllInfected;CA12_5_BinomInitCondit;Trial1.dat';
        fileQ01 = '';%'2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_1_0.5265;MVM_p=0.63_0.03_0.84;Alpha=0.005_0.002_0.017;Alpha_hung=0_1_0;Tau=4_2_20;HungRatio=0.001_0.03_0.21;NumRuns=5000;NumTrials=32;Q0=0.1;CompleteGraph;AllInfected;CA12_5_BinomInitCondit;Trial1.dat';
        fileQ03 = '';%'2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_1_0.5265;MVM_p=0.63_0.03_0.84;Alpha=0.002_0.001_0.007;Alpha_hung=0_1_0;Tau=1_1_8;HungRatio=0.001_0.03_0.21;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_5_BinomInitCondit;Trial1.dat';
        %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_1_0.5265;MVM_p=0.63_0.03_0.84;Alpha=0.001_0.001_0.007;Alpha_hung=0_1_0;Tau=0.5_0.5_5;HungRatio=0.001_0.03_0.21;NumRuns=5000;NumTrials=32;Q0=0.5;CompleteGraph;AllInfected;CA12_5_BinomInitCondit;Trial1.dat';

        %Dts
        %fileQ05 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_1_0.5265;MVM_p=0.63_0.03_0.84;Alpha=0.008_0.004_0.028;Alpha_hung=0_1_0;Tau=1_1_8;HungRatio=0.001_0.03_0.21;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_5_BinomInitCondit.dat';
        %fileQ01 = '2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5265_1_0.5265;MVM_p=0.63_0.03_0.84;Alpha=0.0005_0.0003_0.0018;Alpha_hung=0_1_0;Tau=1_1_8;HungRatio=0.001_0.03_0.21;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;CA12_5_BinomInitCondit.dat';
		
    end
    
    if sum(reshape(data,numel(data),1)) > mindat
        b = 0:length(data(1,:));
        h = data;
        [R,Rp] = bCompareQ0s(h, b,fileQ01,fileQ03,fileQ05);
    	all_R(count,:) = R;
    	all_Rp(count,:) = Rp;
    end
    save(strcat(directory,'CA12VoteTime_New_Dt=15,60,240-',num2str(i),'.mat'),'all_R','all_Rp');

end
%save(strcat(directory,'CA12VoteTime_New_Q0=0.1,0.3,1.mat'),'all_R','all_Rp');

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));

%save(strcat(directory,'Influence_New_Q0=0.1,0.3,1.mat'),'all_R','all_Rp');
%}
end
