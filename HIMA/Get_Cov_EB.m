function [cov_mean,cov_mode]=Get_Cov_EB(W0, n_terms)
    %%%% This function estimates the posterior mode of covaraince matrix
    %%%% from the multivaraite normal distrition using Empirical Bayesian
    %%%% method proposed in [Champion, 2003].
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs:
    %%%%% W0: an nxp matrix containing values of p variables across n
    %%%%% subjects
    %%%%% n_terms, an arguments in the hypergeometirc function below, which
    %%%%% defineds the number of expansion terms needed.
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outputs:    
    %%%%%  cov_mean:   estimated posterior mean 
    %%%%%  cov_mode    estimated posterior mode      


    cov_W0=cov(W0, 'partialrows'); 
    corr_W0=corrcoef(W0,'Rows','pairwise');
    n=size(W0,1); % # of subjects
    p=size(W0,2); % # of variables
 

    % Step1: scalar variances along the diagonals of Z (covariance matrix) 
    z_cov=diag(diag(cov_W0));

    % Step2: estimate off-diagonal elements of Z
    rho_sum=0;
    count=0;
    for i=1:size(W0,2)
        for j=1:size(W0,2)
            if i~=j
                %hypergeometric2f1(a,b,c,x,n_terms), where n_terms is # of expansion terms
                hyper_f1=hypergeometric2f1(0.5,0.5,(n-1)/2,1-corr_W0(i,j)^2 ,n_terms);
                rho_sum=rho_sum+corr_W0(i,j)*hyper_f1; 
                count=count+1;
            end 
        end 
    end 
    rho_mean=rho_sum/count;

    for i=1:size(W0,2)
        for j=1:size(W0,2)
            if i~=j
                z_cov(i,j)=rho_mean*sqrt(  z_cov(i,i)*z_cov(j,j)  );
            end 
        end 
    end 
    
    % Step3: estimate lambda
    %lam=1; %temp
    k2_num=0; %numerator for k^2
    k2_den=0; %denominator for k^2
    
    %correlation matrix of Z_cov
    z_corr=rho_mean;

    for i=1:size(W0,2)
        for j=1:size(W0,2)
            if i~=j
                hyper_f1=hypergeometric2f1(0.5,0.5,(n-1)/2,1-corr_W0(i,j)^2 ,n_terms);
                alpha=corr_W0(i,j)*hyper_f1;

                hyper_f2=hypergeometric2f1(1,1,(n+1)/2,1-corr_W0(i,j)^2 ,n_terms);
                beta=1-(n-2)*(1-corr_W0(i,j)^2)/(n-1)*hyper_f2;
    
                k2_num=k2_num+(beta-2*alpha*z_corr+z_corr^2);
                k2_den=k2_den+ (1-z_corr^2)^2;

            end 
        end 
    end 

    k2=k2_num/k2_den;
    lam=1/k2-3;

    
    % Step4: bayes esitmator   
    cov_mean=(lam.*z_cov+cov_W0)/(lam+n);

    cov_mode=(lam.*z_cov+cov_W0)/(lam+n+2*p+2);

%     % Step5: adjustment(rescale cov_bayes_W0):
%     for i=1:size(W0,2)
%         for j=1:size(W0,2)
%             if i~=j
%                 cov_bayes_W0(i,j)=cov_bayes_W0*sqrt(cov_W0(i,i)*cov_W0(j,j)/(z(i,i)*z(j,j)));
%             end 
%         end 
%     end 

    % Check if covariance matrix SPD
    d=eig(cov_mean); %eigenvalues
    if all(d>0)==0 %if non-PD
        cov_mean=nearestSPD(cov_mean);
    end 

    d=eig(cov_mode); %eigenvalues
    if all(d>0)==0 %if non-PD
        cov_mode=nearestSPD(cov_mode);
    end

end 
