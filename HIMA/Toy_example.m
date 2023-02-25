% This is a top example to implement HIMA to impute a synthetic data. 


%% Kernel smooth on subjects: 
% align their smoothed probability density estimate so that the subject-level 
% random effect is taken into account.

%W0: input matrix containing values of p variables across n subjects
%import W0

W0_mean = mean(W0.', 'omitnan'); %1*58
total_mean = mean(W0_mean);
W0 = W0 - (W0_mean - total_mean).';


%% Imputation using HIMA
inv_lambda = 0.1; %Tuning parameters for the diagonals of covariance
num_imp = 15;    % the number of imputed datasets (M) 
num_iter = 20;   % the number of interations (T)

tic
[imp_res, mu_all,cov_all] =Get_Imp_EB2(W0,inv_lambda,num_imp, num_iter);
    %%%%%  imp_res:  imputation results, which is a 1xM cell. Each cell 
    %%%%% (i.e., imp_res{1,m}) contains an imputed nxp matrix.

    %%%%%  mu_all: a 1xM cell. Each cell contains a pxT matrix, recording 
    %%%%%  the updated px1 mu vector over T iterations. 
    %%%%%  mu_all can be used to plot traceplot and examine convegence.

    %%%%%  cov_all: a MxT cell. Each cell contains a pxp covariance matrix
    %%%%%  computed at the t-th iteration in the m-th imputed dataset.
    %%%%%  cov_all can be used to plot traceplot and examine convegence.
timeElapsed = toc %running time in seconds


%% Scaled back to imputed results
imp_res = cellfun(@(x) x + (W0_mean - total_mean).', imp_res_scaled, 'UniformOutput', false);

%% Traceplot of mean and covaraince 

voxel_j=120; %select a voxel (e.g. the 120-th voxel) to plot its traceplot
x=1:num_iter; 

figure;
set(gca,'FontSize',18)
set(gca,'YTick',64.76:0.01:64.8)
for m=0:4
    mu_m=mu_all(:,m*iter+1:(m+1)*iter);
    y=mu_m(voxel_j,:);
    
    hold on
    plot(x,y,'-','LineWidth',1.5)   
end 
% xlabel('Iteration numbers','FontSize',16,'FontWeight','bold')
% ylabel('Mean of Voxel #120 ','FontSize',16,'FontWeight','bold')

%title('Our method','fontweight','bold','fontsize',14)
% title('Trace plot of mean in MVN for voxel 120','fontweight','bold','fontsize',14)
hold off


%% Inspection 
% Inspect the distributions of the original values and the imputed values
% Use voxel #120 and #305 as an example

check_missingrate = mean(isnan(W0))'; 
test=sort(check_missingrate,"descend");
high_idx=find(check_missingrate<0.2&check_missingrate>0.17);
check_missingrate(high_idx);
[high_idx check_missingrate(high_idx)]  %for easier visualization

%%%%%%%%%%% Scatter plots %%%%%%%%%%
figure
subplot(1,2,1)
j=120
mis=imp_res(:,j);obs=W0(:,j);
DataArray=nan(size(obs,1),2);
DataArray(1:length(mis),1)=mis;DataArray(1:length(obs),2)=obs;
UnivarScatter(DataArray,'Label',{'imputed','original'});
title("Voxel "+j,'fontweight','bold','fontsize',14)

j=305
mis=imp_res(:,j);obs=W0(:,j);
DataArray=nan(size(obs,1),2);
DataArray(1:length(mis),1)=mis;DataArray(1:length(obs),2)=obs;
UnivarScatter(DataArray,'Label',{'imputed','original'});
title("Voxel "+j,'fontweight','bold','fontsize',14)

sgtitle('Right hippocampus','fontweight','bold','fontsize',16) 


%%%%%%%%%% density plots %%%%%%%%%%
figure;
subplot(1,2,1) 
j=120; %voxel
obs=X_w_NA_saveMR(:,j);

hold on
m=randsample(num_imp,50);  %unique 

for i=1:length(m)
    [f,xi] = ksdensity(imp_res_saveMR{1,m(i)}(:,j));
    plot(xi,f,'LineWidth',3);
end 
[f1,xi1] = ksdensity(obs);
plot(xi1,f1,'LineWidth',3);

title("Voxel "+j,'fontweight','bold','fontsize',14)
legend({'imputed','original'},'Location','best')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2) 
j=305; %voxel
obs=X_w_NA_saveMR(:,j);

hold on
% num_imp = 100;    %% number of imputation datasets
m=randsample(num_imp,50);  %unique 

for i=1:length(m)
    [f,xi] = ksdensity(imp_res_saveMR{1,m(i)}(:,j));
    plot(xi,f,'LineWidth',3);
end 

[f1,xi1] = ksdensity(obs);
plot(xi1,f1,'LineWidth',3);
%xlim([0, 200])
%ylim([0, 0.0122])
title("Voxel "+j,'fontweight','bold','fontsize',14)
legend({'imputed','original'},'Location','best')
hold off
sgtitle('Right hippocampus','fontweight','bold','fontsize',16) 







