(* ::Package:: *)

function [sumLs] = all_data _parse _vote _NewInfluence _TTNull()
directory = '/Users/keithburghardt/Google Drive/DelibWork/DelibData/';
load(strcat(directory,'VoteBinnedData_AllStates _CrimCivil.mat'));

data = eval('OR6Civil');
b = (0:length(data(1,:)))./60;
h = data;
N = numel(h(:,1))-1;
LambdaRangeGuilty=0.1:0.1:2;LambdaRangeInnocent=0.1:0.1:2;
[y, t] = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h, b);
[param_est, L,sumL] = bTTfit(data,b,y,LambdaRangeGuilty,LambdaRangeInnocent);
OR6CivilsumL=sumL;


data = eval('OR6Criminal');
b = (0:length(data(1,:)))./60;
h = data;
N = numel(h(:,1))-1;
LambdaRangeGuilty=0.1:0.1:2;LambdaRangeInnocent=0.1:0.1:2;
[y, t] = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h, b);
[param_est, L,sumL] = bTTfit(data,b,y,LambdaRangeGuilty,LambdaRangeInnocent);
OR6CriminalsumL=sumL;


data = eval('OR12Civil');
b = (0:length(data(1,:)))./60;
h = data;
N = numel(h(:,1))-1;
LambdaRangeGuilty=0.1:0.1:2;LambdaRangeInnocent=0.1:0.1:2;
[y, t] = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h, b);
[param_est, L,sumL] = bTTfit(data,b,y,LambdaRangeGuilty,LambdaRangeInnocent);
OR12CivilsumL=sumL;


data = eval('OR12Criminal');
b = (0:length(data(1,:)))./60;
h = data;
N = numel(h(:,1))-1;
LambdaRangeGuilty=0.1:0.1:2;LambdaRangeInnocent=0.1:0.1:2;
[y, t] = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h, b);
[param_est, L,sumL] = bTTfit(data,b,y,LambdaRangeGuilty,LambdaRangeInnocent);
OR12CriminalsumL=sumL;

sumLs = [OR6CivilsumL OR6CriminalsumL OR12CivilsumL OR12CriminalsumL];
end
