function [CDFOR6,CDFOR12,CDFCA6,CDFCA8,CDFCA12] = all_data_parse_vote_TTExPlots() 

directory = '/export/data/ccbdata/keith/DelibWork/DelibData/'; 
load(strcat(directory,'NewBin_VoteBinnedData.mat'));
%load(strcat(directory,'OR6VoteTimeDistModel_NullInfluence.mat'));% CHANGE
load(strcat(directory,'Influence_Null_5_17_16_1minT_NewCA.mat'));% CHANGE
mindat = 23;
LambdaRangeGuilty = 0.1:0.05:3; 
LambdaRangeInnocent = 0.1:0.05:3;

count = 0;

low = 1;
high = 2;


%OR 6 
for i=1:1
    count = count + 1;

    data = eval(strcat('VoteBinnedOR6data',num2str(i)));
    
    b = (0:length(data(1,:)))./60;
    h = data;
    N = numel(h(:,1))-1;
    TTCDF = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h,b);
    CDFOR6 = GetTTCDF(h,b,TTCDF,LambdaRangeGuilty,LambdaRangeInnocent);
    Params = reshape(all_param_quantiles(count,8:9,2:3),2,2);
    CDFOR6_Low  =  FindCDF_TT(h,TTCDF,Params(1,low), Params(2,low),LambdaRangeGuilty,LambdaRangeInnocent);
    CDFOR6_High =  FindCDF_TT(h,TTCDF,Params(1,high),Params(2,high),LambdaRangeGuilty,LambdaRangeInnocent);

    save(strcat(directory,'OR6VoteTimeTTBestFitModel_1minT.mat'),'CDFOR6','CDFOR6_Low','CDFOR6_High');
end

%OR 12 
for i=1:1
    count = count + 1;
    
    data = eval(strcat('VoteBinnedOR12data',num2str(i)));
    
    b = (0:length(data(1,:)))./60;
    h = data;
    TTCDF = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h,b);
    CDFOR12 = GetTTCDF(h,b,TTCDF,LambdaRangeGuilty,LambdaRangeInnocent);
    Params = reshape(all_param_quantiles(count,8:9,2:3),2,2)
    CDFOR12_Low  =  FindCDF_TT(h,TTCDF,Params(1,low), Params(2,low),LambdaRangeGuilty,LambdaRangeInnocent);
    CDFOR12_High =  FindCDF_TT(h,TTCDF,Params(1,high),Params(2,high),LambdaRangeGuilty,LambdaRangeInnocent);

    save(strcat(directory,'OR12VoteTimeTTBestFitModel_1minT.mat'),'CDFOR12','CDFOR12_Low','CDFOR12_High'); 
end

%CA 6 
    n = 6;
    data = zeros(n+1,numel(VoteBinnedCA6data(1,1,:)));
    for i=1:9
        count = count + 1;
    	data = data + reshape(VoteBinnedCA6data(i,:,:),n+1,numel(VoteBinnedCA6data(i,1,:)));
    end
    if sum(reshape(data,numel(data),1)) > mindat
        b = 0:length(data(1,:));
        %No12Hr = b(1:end-1)==12;
        h = data;
        %h(:,No12Hr) = 0;
        TTCDF = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h,b);
        CDFCA6 = GetTTCDF(h,b,TTCDF,LambdaRangeGuilty,LambdaRangeInnocent);
        Params = reshape(all_param_quantiles(count,8:9,2:3),2,2);
        CDFCA6_Low  =  FindCDF_TT(h,TTCDF,Params(1,low), Params(2,low),LambdaRangeGuilty,LambdaRangeInnocent);
        CDFCA6_High =  FindCDF_TT(h,TTCDF,Params(1,high),Params(2,high),LambdaRangeGuilty,LambdaRangeInnocent);
    end
    save(strcat(directory,'CA6VoteTimeTTBestFitModel_1minT.mat'),'CDFCA6','CDFCA6_Low','CDFCA6_High'); 


%CA8 
n = 8; 
data = reshape(VoteBinnedCA8data(1,:,:),n+1,numel(VoteBinnedCA8data(1,1,:)));
b = 0:length(data(1,:));
CDFCA8 = zeros(9,n+1,numel(b));
CDFCA8_Low = zeros(9,n+1,numel(b));
CDFCA8_High = zeros(9,n+1,numel(b));

for i=1:9
    count = count + 1;
    data = reshape(VoteBinnedCA8data(i,:,:),n+1,numel(VoteBinnedCA8data(i,1,:)));
    
    if sum(reshape(data,numel(data),1)) > mindat
        b = (0:length(data(1,:)));
        %No12Hr = b(1:end-1)==12;
        h = data;
        %h(:,No12Hr) = 0;
	TTCDF = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h,b);
        CDFCA8(i,:,:) = GetTTCDF(h,b,TTCDF,LambdaRangeGuilty,LambdaRangeInnocent);
        Params = reshape(all_param_quantiles(count,8:9,2:3),2,2);
        %CDFCA8_Low(i,:,:)  =  FindCDF_TT(h,TTCDF,Params(1,low), Params(2,low), LambdaRangeGuilty,LambdaRangeInnocent);
        %CDFCA8_High(i,:,:) =  FindCDF_TT(h,TTCDF,Params(1,high),Params(2,high),LambdaRangeGuilty,LambdaRangeInnocent);
    end
end 
save(strcat(directory,'CA8VoteTimeTTBestFitModel_1minT.mat'),'CDFCA8','CDFCA8_Low','CDFCA8_High'); 

%CA12 
n=12; 
data = reshape(VoteBinnedCA12data(1,:,:),n+1,numel(VoteBinnedCA12data(1,1,:)));
b = 0:length(data(1,:));
CDFCA12 = zeros(9,n+1,numel(b)); 
CDFCA12_Low = zeros(9,n+1,numel(b)); 
CDFCA12_High = zeros(9,n+1,numel(b)); 

for i=1:9
    count = count + 1;
    data = reshape(VoteBinnedCA12data(i,:,:),n+1,numel(VoteBinnedCA12data(i,1,:)));
    
    if sum(reshape(data,numel(data),1)) > mindat
        b = 0:length(data(1,:));
        %No12Hr = b(1:end-1)==12;
        h = data;
        %h(:,No12Hr) = 0;
	TTCDF = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h,b);
        CDFCA12(i,:,:) = GetTTCDF(h,b,TTCDF,LambdaRangeGuilty,LambdaRangeInnocent);
        Params = reshape(all_param_quantiles(count,8:9,2:3),2,2);
        %CDFCA12_Low(i,:,:)  =  FindCDF_TT(h,TTCDF,Params(1,low), Params(2,low),LambdaRangeGuilty,LambdaRangeInnocent);
        %CDFCA12_High(i,:,:) =  FindCDF_TT(h,TTCDF,Params(1,high),Params(2,high),LambdaRangeGuilty,LambdaRangeInnocent);
    end
end 
save(strcat(directory,'CA12VoteTimeTTBestFitModel_1minT.mat'),'CDFCA12','CDFCA12_Low','CDFCA12_High'); 

save(strcat(directory,'TTBestFitModels_1minT.mat'),'CDFOR6','CDFOR12','CDFCA6','CDFCA8','CDFCA12','CDFOR6_Low','CDFOR12_Low','CDFCA6_Low','CDFCA8_Low','CDFCA12_Low','CDFOR6_High','CDFOR12_High','CDFCA6_High','CDFCA8_High','CDFCA12_High'); 

end
