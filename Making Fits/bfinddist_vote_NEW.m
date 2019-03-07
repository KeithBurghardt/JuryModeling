function [p, param_est,param_errors,param_quantiles, N_tail,Bmin,Bmin_quantiles] = bfinddist_vote_NEW(h, Tboundaries,dt,file,InfluenceCDF,reps)
    % What this function does:
    %   - Find estimates of values, error bars, & bmin (if exponential) 
    %   - Find p-value of each function
    %   - If p-value > 0.1, compare noramalized likelihood ratios
    %

    %Total amount of data
    n = sum(reshape(h,numel(h),1));

    % number of jurors
    N = length(h(:,1)) - 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   Majority Voter Model - Stubbornness    %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Influence_param_est = zeros(1,6)-1;
    Influence_param_errors = zeros(1,6)-1;
    Influence_param_quantiles = zeros(6,4)-1;
    p_Influence = -1;
   
    %[Params, BestL] = bInfluencefit(h,Tboundaries,InfluenceCDF,file);   
    %disp('BEST FIT:'); 
    %disp(Params);
    %disp('BEST L 10^-11:'); 
    %disp(BestL);

    % find estimates, error bars
    [Influence_param_errors, eof] = bInfluencevar(h,Tboundaries,InfluenceCDF,file,'reps',reps);

    % estimated lambda for evl
    Influence_param_est  = [  mean(eof(:,1)) ...
                mean(eof(:,2)) ...
                mean(eof(:,3)) ...
                mean(eof(:,4)) ...
                mean(eof(:,5)) ...
                mean(eof(:,6))
                ];

    % quantiles of lambda
    Influence_param_quantiles = [quantile(eof(:,1),[0.025 0.05 0.95 0.975]);
                        quantile(eof(:,2),[0.025 0.05 0.95 0.975]);
                        quantile(eof(:,3),[0.025 0.05 0.95 0.975]);
                        quantile(eof(:,4),[0.025 0.05 0.95 0.975]);
                        quantile(eof(:,5),[0.025 0.05 0.95 0.975]);
                        quantile(eof(:,6),[0.025 0.05 0.95 0.975])
                        ];
    % p-value
    p_Influence = bInfluencepva(h, Tboundaries,InfluenceCDF,file,'reps',reps);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%                 TT_Null                  %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    LambdaRangeGuilty = 0.1:0.05:3;
    LambdaRangeInnocent = 0.1:0.05:3;
    TT_param_est = zeros(1,2)-1;
    TT_param_errors = zeros(1,2)-1;
    TT_param_quantiles = zeros(2,4)-1;
    p_TT = -1;
    
    %{    
    TTCDF = TTcdf_paramspace(LambdaRangeGuilty,LambdaRangeInnocent,h, Tboundaries);

    % find estimates, error bars
    [TT_param_errors, eof] = bTTvar(h,Tboundaries,TTCDF,LambdaRangeGuilty,LambdaRangeInnocent, 'reps',reps);

    % estimated lambda for evl
    TT_param_est  = [  mean(eof(:,1)) ...
                    mean(eof(:,2))           
                    ];

    % quantiles of lambda
    TT_param_quantiles = [quantile(eof(:,1),[0.025 0.05 0.95 0.975]);
        quantile(eof(:,2),[0.025 0.05 0.95 0.975])
        ];
    % p-value
    p_TT =  bTTpva(h,Tboundaries,TTCDF,LambdaRangeGuilty,LambdaRangeInnocent,'reps',reps);
    %}  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%         Exponential             %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: split votes, find exponential tail of each final vote
    % Record: final vote, t_min, exponent, error, KS test p-value
    % set values to 0
    N_tail = zeros(2,N+1)-1;
    lambda_est = zeros(1,N+1)-1;
    sigma_lambda = zeros(1,N+1)-1;
    Bmin = zeros(1,N+1);
    Bmin_quantiles = zeros(N+1,4);
    %show 90, 95% confidence intervals
    lambda_quantile = zeros(N+1,4)-1;
    p_exp = zeros(1,N+1)-1;
    %{
    if length(Tboundaries) > 100
        limit = 2; %exponential slope is seen before 2 hours
    else
        limit = 8; %exponential slope is seen before 8 hours
    end
    
    for j = 1:N+1
        if sum(h(j,:))>20
            % find estimates, error bars, bmin

            [sigma_lambda(j), ~, ~, eof] = bexpvar(h(j,:),Tboundaries,'limit',limit,'reps',reps);

            bmin_est = mean(eof(:,2));
            bmin_quantiles = quantile(eof(:,2),[0.025 0.05 0.95 0.975]);
            % find bmin closest to a boundary value
            [~,pos] = min(abs(Tboundaries - bmin_est));
            bmin_est = Tboundaries(pos);
            % amount of data at the tail
            ind = find(Tboundaries>=bmin_est, 1);
            Bmin(j) = bmin_est;
            Bmin_quantiles(j,:) = bmin_quantiles;
            N_tail(:,j) = [sum(h(j,:)); sum(h(j,ind:end))];
            %N_tot = [N N_tail];
            % estimated lambda
            lambda_est(j) = mean(eof(:,3));
            lambda_quantile(j,:) = quantile(eof(:,3),[0.025 0.05 0.95 0.975]);
            fprintf('Exponential:\n lambda = %f +/- %f\n bmin = %f\n\n', lambda_est(j), sigma_lambda(j),bmin_est);
            if N_tail(2,j) > 25   
                
                % p-value, using the above estimated bmin
                p_exp(j) = bexppva(h(j,:), Tboundaries, bmin_est,'limit',limit,'reps',reps);
                fprintf('p = %f\n\n', p_exp(j));
            end
        end
        
    end
    %}    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%        Final Vote Null          %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p_VoteNull = -1;
    
    EmpiricalVote = sum(h,2);
    VoteNullParams = zeros(1,2) - 1;
    p_VoteNull = -1;
    %[p_VoteNull,~,~,VoteNullParams] = NullVotepva(EmpiricalVote,'reps',reps);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%    Comparisons of Distributions   %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ALL p-values
    % 1 x N+4
    p = [p_Influence p_VoteNull p_TT p_exp];

    % ALL parameter estimates. 
    % 4 x N+1
    param_est = [Influence_param_est zeros(1,N+1-6); VoteNullParams zeros(1,N+1-2); TT_param_est zeros(1,N+1-2); lambda_est];

    % ALL parameter errors.
    % 4 x N+1
    param_errors = [Influence_param_errors zeros(1,N+1-6); zeros(1,N+1); TT_param_errors zeros(1,N+1-2);sigma_lambda];

    % ALL parameter quantiles
    % (N + 10) x 4
    param_quantiles = [Influence_param_quantiles; zeros(1,4); TT_param_quantiles; lambda_quantile];


end
