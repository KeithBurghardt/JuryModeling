function [R,p] = bCompareQ0s(h, Tboundaries,fileQ01,fileQ03,fileQ05)
    % What this function does:
    %   - Find estimates of values, error bars, & bmin (if exponential) 
    %   - Find p-value of each function
    %   - If p-value > 0.1, compare noramalized likelihood ratios
    %

    %Total amount of data
    n = sum(reshape(h,numel(h),1));

    % number of jurors
    N = length(h(:,1)) - 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   Majority Voter Model - Stubbornness (Q0 = 0.1)   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = 60;%change to 15
    N = numel(h(:,1)) - 1;

    InfluenceCDF = Influencecdf_paramspace(fileQ01,dt,Tboundaries);
    [~, lMVM_01] = bInfluencefit(h, Tboundaries,InfluenceCDF,fileQ01);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   Majority Voter Model - Stubbornness (Q0 = 0.3)   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = 60;
    N = numel(h(:,1)) - 1;

    InfluenceCDF = Influencecdf_paramspace(fileQ03,dt,Tboundaries);
    [~, lMVM_03] = bInfluencefit(h, Tboundaries,InfluenceCDF,fileQ03);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   Majority Voter Model - Stubbornness (Q0 = 0.5)   %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = 60;% change to 60*4
    N = numel(h(:,1)) - 1;
    InfluenceCDF = Influencecdf_paramspace(fileQ05,dt,Tboundaries);
    [~, lMVM_05] = bInfluencefit(h, Tboundaries,InfluenceCDF,fileQ05);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%    Comparisons of Distributions   %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [R_MVM_01_03,p_MVM_01_03] = blrtest_l1l2(lMVM_01,lMVM_03,0);
    [R_MVM_01_05,p_MVM_01_05] = blrtest_l1l2(lMVM_01,lMVM_05,0);
    [R_MVM_03_05,p_MVM_03_05] = blrtest_l1l2(lMVM_03,lMVM_05,0);


    R = [R_MVM_01_03 R_MVM_01_05 R_MVM_03_05];
    p = [p_MVM_01_03 p_MVM_01_05 p_MVM_03_05];

end
