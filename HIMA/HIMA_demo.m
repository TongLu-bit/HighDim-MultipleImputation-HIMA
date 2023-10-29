%% HIMA: A Multiple Imputation Method for High-Dimensional Neuroimaging Data
% HIMA develops a new computational strategy for sampling large covariance matrices 
% based on a robustly estimated posterior mode, 
% significantly enhancing computational efficiency and numerical stability.  

%% Import data
% W0: the input matrix containing values of p variables across n subjects
load('data.mat');
W0=data;
%% Kernel smooth on subjects: 
% Align smoothed probability density estimate so that the subject-level 
% random effect is taken into account.

W0_mean = mean(W0.', 'omitnan'); %1*58
total_mean = mean(W0_mean);
W0 = W0 - (W0_mean - total_mean).';

%% Imputation using HIMA
% The main result that the function "HIMA" returns is called
% "imp_res_scaled", which is a 1xM cell. Each cell (i.e., imp_res{1,m}) contains an imputed nxp matrix.
inv_lambda = 0.1; % Tuning parameters for the diagonals of covariance
num_imp = 5;     % the number of imputed datasets (M) 
num_iter = 30;    % the number of interations (T)
if_store_mu=1;   % 0 if don't want to store the updated mu over iterations
if_store_cov=0;  % 0 if don't want to store the updated cov over iterations

tic
[imp_res_scaled, mu_all,cov_all] =HIMA(W0,inv_lambda,num_imp, num_iter,if_store_mu, if_store_cov);
    %%%%%  imp_res_scaled:  imputation results, which is a 1xM cell. Each cell 
    %%%%% (i.e., imp_res{1,m}) contains an imputed nxp matrix.

    %%%%%  mu_all: a 1xM cell. Each cell contains a pxT matrix, recording 
    %%%%%  the updated px1 mu vector over T iterations. 
    %%%%%  mu_all can be used to plot traceplot and examine convegence.

    %%%%%  cov_all: a MxT cell. Each cell contains a pxp covariance matrix
    %%%%%  computed at the t-th iteration in the m-th imputed dataset.
    %%%%%  cov_all can be used to plot traceplot and examine convegence.
timeElapsed = toc %running time in seconds
display("in seconds")

%% Scaled back to imputed results
% Reverse back to the scale before performing kernel smoothing
imp_res = cellfun(@(x) x + (W0_mean - total_mean).', imp_res_scaled, 'UniformOutput', false);

%% Traceplot of mean 
voxel_j=50; %select a voxel (e.g. the 120-th voxel) to plot its traceplot
x=1:num_iter; 

set(figure, 'Position', [100, 100, 800, 400]);
set(gca,'FontSize',18)
set(gca,'YTick')
for m=1:num_imp % go through each imputed dataset 
    mu_dataset=mu_all{1,m};
    y=mu_dataset(voxel_j,:);
    
    hold on
    plot(x,y,'-','LineWidth',1.5)   
end 
xlabel('Iteration numbers','FontSize',16,'FontWeight','bold')
ylabel('Mean of Voxel #120 ','FontSize',16,'FontWeight','bold')
title('Trace plot of mean in MVN for a random selected voxel','fontweight','bold','fontsize',14)
hold off

