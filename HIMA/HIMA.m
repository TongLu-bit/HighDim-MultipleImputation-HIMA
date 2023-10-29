function [imp_res, mu_all,cov_all] =HIMA(W0,inv_lambda,num_imp, num_iter,if_store_mu, if_store_cov)
    %%%% This function implements the HIMA mode to impute missingness in
    %%%% high-dementional data W0, under the assumption of MAR and the
    %%%% architecture of MCMC.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs:
    %%%%% W0: an nxp matrix containing values of p variables across n
    %%%%% subjects.
    %%%%% inv_lambda: tuning parameters for the diagonals of covariance
    %%%%% matrix.
    %%%%% num_imp: the number of imputed datasets (denoted by M) .
    %%%%% num_iter: the number of interations (denoted by T).
    %%%%% if_store_mu: 0 if don't want to store the updated px1 mu matrix 
    %%%%% at each iteration (for convergence check purpose); 1 otherwise.
    %%%%% if_store_cov: 0 if don't want to store the updated big pxp cov matrix 
    %%%%% at each iteration (for convergence check purpose); 1 otherwise.
    %%%%% Generaly, when p is large, it is recommended to set if_store_cov
    %%%%% to 0 as it will occupy memory and largely slow down running time.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outputs:    
    %%%%%  imp_res:  imputation results, which is a 1xM cell. Each cell 
    %%%%% (i.e., imp_res{1,m}) contains an imputed nxp matrix.
    %%%%%  mu_all: a 1xM cell if the argument if_store_mu=1. Each cell contains a pxT matrix, recording 
    %%%%%  the updated px1 mu vector over T iterations. 
    %%%%%  mu_all can be used to plot traceplot and examine convegence.
    %%%%%  cov_all: a MxT cell if the argument if_store_cov=1. Each cell contains a pxp covariance matrix
    %%%%%  computed at the t-th iteration in the m-th imputed dataset.
    %%%%%  cov_all can be used to plot traceplot and examine convegence.

    numRow_W0 = size(W0, 1);
    numCol_W0 = size(W0, 2);

    %%%%%%%% Initialization %%%%%%%%%%%%%% 
    % Initial imputation (using mean inputation)
    W=W0;
    mean_vector=mean(W0,'omitnan');
    for j=1:numCol_W0
        mis_idx=find(isnan(W0(:,j))==1);    
        W(mis_idx,j)=mean_vector(j);
    end 

    % Initialize mean and covariance matrix:
    mu_W=mean(W,1)';
    cov_W=cov(W);
    d=eig(cov_W); %eigenvalues
    if all(d>0)==0 %if non-PD
        cov_W=nearestSPD(cov_W);
    end

    %%%%%%%%% MCMC process %%%%%%%%%%%%%%
    
    imp_res = cell(1, num_imp); % store imputation result
    mu_all=cell(1, num_imp);  %store updated mu matrix over iterations
    cov_all=cell(num_imp, num_iter); %store updated cov matrix over iterations
    
    for m = 1:num_imp  
         
        text1 = ['Current imputed data set: m=',num2str(m)];
        disp(text1)

        W0_imp = W;  % Reset to default W for each imputation dataset. 

        for t = 1:num_iter

            if mod(t,10)==0
                text2 = ['Current iteration: t=',num2str(t)];
                disp(text2)
            end

            for s = 1:numRow_W0
                     
                mis_idx = find(isnan(W0(s,:)) == 1); 

                if ~isempty(mis_idx) %skip if the current subject does not contain missingnes
                    
                    obs_idx = setdiff(1:numCol_W0, mis_idx);

                    %%%%%%%%%%%%%% Imputation I-step: generate Y %%%%%%%%%%%%%%

                    % Pre-compute the intermediate term, Sigma_12*sigma_22^-1, to speed up computaiton
                    %inter_term=cov_bayes_W0(mis_idx, obs_idx)/  ( cov_bayes_W0(obs_idx, obs_idx) + inv_lambda * eye(length(obs_idx))   );
                    %inter_term=cov_W(mis_idx, obs_idx)/ cov_W(obs_idx, obs_idx);
                    inter_term=cov_W(mis_idx, obs_idx)/ ( cov_W(obs_idx, obs_idx)+ inv_lambda * eye(length(obs_idx)) );

                    % calculate mu_(mis|obs):
                    mu_mis_given_obs=mu_W(mis_idx)+inter_term*(W0_imp(s, obs_idx)' - mu_W(obs_idx));

                    % calculate sigma_(mis|obs):
                    cov_mis_given_obs=cov_W(mis_idx, mis_idx) - inter_term * cov_W(obs_idx, mis_idx);
                    cov_mis_given_obs = (cov_mis_given_obs +cov_mis_given_obs.')/2;   %% To avoid non-symmetric issue

                    % impute Y_mis
                    na_imp = mvnrnd(mu_mis_given_obs, cov_mis_given_obs, 500); 
                    na_imp_mean=mean(na_imp); 

                    % Replace Y=(Y_obs, Y_mis^[t])
                    W0_imp(s, mis_idx) = na_imp_mean;

                    %%%%%%%%%%%%%% Posterior P-step  %%%%%%%%%%%%%%
                    % Draw mu
                    pos_mu=mean(W0_imp, 1)';
                    pos_sigma=cov_W/length(pos_mu);

                    mu_W = mvnrnd(pos_mu, pos_sigma, 1);
                    mu_W=mu_W';
                    %mu_W=mean(mu_W)';
                    
                    % Approximate posterior mode of covariance
                    [~,cov_W]=Get_Cov_EB(W0_imp, 50);                  
                    
                end
            end

            if if_store_mu==1
                mu_all{1,m}(:,t)=mu_W;   
            end 
            if if_store_cov==1
                cov_all{m,t}=cov_W; 
            end

        end

        imp_res{1, m} = W0_imp;
        if if_store_mu==0
            disp("Updated normal mean is chosen not to be recorded over iterations");
        end 
        
        if if_store_cov==0
            disp("Updated normal covariance is chosen not to be recorded over iterations");
        end 
     
    end
end

   

