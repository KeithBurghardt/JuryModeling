function [all_Rp,all_R] = all_data_parse_vote_NewInfluence_Compare()
%directory = '/export/data/ccbdata/keith/DelibWork/DelibData/';
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
    fileDt15='';
    fileDt60='';
    fileDt240='';
    [R,Rp] = bCompareDts(h, b,fileDt15,fileDt60,fileDt240);
    save(strcat(directory,'OR6Civil_Comparison_15_60_240s.mat'),'R','Rp');
    all_R(count,:) = R;
    all_Rp(count,:) = Rp;
end

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));

%OR 6 Criminal
for i=1:1
    count = count + 1;

    %data = eval(strcat('VoteBinnedOR6data',num2str(i)));
    data = eval(strcat('OR6Criminal'));

    b = (0:length(data(1,:)))./60;
    h = data;
    fileDt15='';
    fileDt60='';
    fileDt240='';
    [R,Rp] = bCompareDts(h, b,fileDt15,fileDt60,fileDt240);
    save(strcat(directory,'OR6Criminal_Comparison_15_60_240s.mat'),'R','Rp');
    all_R(count,:) = R;
    all_Rp(count,:) = Rp;
end

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));

%OR 12 Civil
for i=1:1
    count = count + 1;
    
    %data = eval(strcat('VoteBinnedOR12data',num2str(i)));
    data = eval(strcat('OR12Civil'));
    
    b = (0:length(data(1,:)))./60;
    h = data;
        fileDt15='';
    fileDt60='';
    fileDt240='';
    [R,Rp] = bCompareDts(h, b,fileDt15,fileDt60,fileDt240);
    save(strcat(directory,'OR12Civil_Comparison_15_60_240s.mat'),'R','Rp');
    all_R(count,:) = R;
    all_Rp(count,:) = Rp;
end

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));

%OR 12 Criminal
for i=1:1
    count = count + 1;
    
    %data = eval(strcat('VoteBinnedOR12data',num2str(i)));
    data = eval(strcat('OR12Criminal'));
    
    b = (0:length(data(1,:)))./60;
    h = data;
        fileDt15='';
    fileDt60='';
    fileDt240='';
    [R,Rp] = bCompareDts(h, b,fileDt15,fileDt60,fileDt240);
    save(strcat(directory,'OR12Criminal_Comparison_15_60_240s.mat'),'R','Rp');
    all_R(count,:) = R;
    all_Rp(count,:) = Rp;
end

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));
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
    fileDt15='';
    fileDt60='';
    fileDt240='';
    [R,Rp] = bCompareDts(h, b,fileDt15,fileDt60,fileDt240);
    all_R(count,:) = R;
    all_Rp(count,:) = Rp;
%end
save(strcat(directory,'CA6VoteTimeDistModel_NullInfluence_15_60_240s_NewCA.mat'),'all_R','all_Rp');

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));
%}

%CA8
for i=1:5
    count = count + 1;
    n = 8;
    data = reshape(VoteBinnedCA8data(i,:,:),n+1,numel(VoteBinnedCA8data(i,1,:)));
    
    if sum(reshape(data,numel(data),1)) > mindat
        b = (0:length(data(1,:)));
        h = data;
        fileDt15='';
        fileDt60='';
        fileDt240='';
        [R,Rp] = bCompareDts(h, b,fileDt15,fileDt60,fileDt240);
    	all_R(count,:) = R;
    	all_Rp(count,:) = Rp;
    end
end
save(strcat(directory,'CA8VoteTimeDistModel_NullInfluence_15_60_240s_NewCA.mat'),'all_R','all_Rp');

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));

%CA12
for i=1:7
    count = count + 1;
    n = 12;
    data = reshape(VoteBinnedCA12data(i,:,:),n+1,numel(VoteBinnedCA12data(i,1,:)));
    
    if sum(reshape(data,numel(data),1)) > mindat
        b = 0:length(data(1,:));
        h = data;
        fileDt15='';
        fileDt60='';
        fileDt240='';
        [R,Rp] = bCompareDts(h, b,fileDt15,fileDt60,fileDt240);
    	all_R(count,:) = R;
    	all_Rp(count,:) = Rp;
    end
end
save(strcat(directory,'CA12VoteTimeDistModel_NullInfluence_15_60_240s_NewCA.mat'),'all_R','all_Rp');

disp(all_R(1:count,:));
disp(all_Rp(1:count,:));

save(strcat(directory,'Influence_Null_15_60_240s_NewCA.mat'),'all_R','all_Rp');

end
