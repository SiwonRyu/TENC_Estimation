clear all
addpath('mfiles');

baseroot = 'C:\rsw\Replication_CP_2021\Empirical Illustration\';
Dyadic_data_raw = importdata([baseroot,'Dyadic_data.xlsx']);
Ind_data_raw = importdata([baseroot,'Ind_data.xlsx']);

Data_dyadic     = Dyadic_data_raw.data;
Data_individual = Ind_data_raw.data;
dyadic_index    = Data_dyadic(:,2)~=Data_dyadic(:,3);
Data_dyadic     = Data_dyadic(dyadic_index,:);

% Load from dyadic data
idx_d     = Data_dyadic(:,1:3);
A0      = Data_dyadic(:,4);
A1      = Data_dyadic(:,5);
DA      = A1-A0;
Di      = Data_dyadic(:,6);
Dj      = Data_dyadic(:,7);

% Load from individual data
idx_i = Data_individual(:,[1,2]);
Y0      = Data_individual(:,3);
Y1      = Data_individual(:,4);
DY      = Y1-Y0;
D       = Data_individual(:,5);
%%
clc
% Sorting for indivudal data is required to correctly construct Z, X matrix
Data_individual = sortrows(Data_individual, [1,2]); % Sort by g -> i

% Random Experiment
N = 47;
Data_dyadic_RE      = [idx_d, Di Dj A1];
Data_individual_RE  = [idx_i, D, Y1];

disp('Estimation Under Randomized Experimental Setting');disp(' ');
[est_zeta1, est_beta_RE, est_pi_RE, est_beta_L_RE, est_beta_S] ...
    = est_RE(Data_dyadic_RE,Data_individual_RE,N);
est_zeta1
est_beta_RE
est_beta_L_RE
est_pi_RE

disp(' ');disp(' ');
disp('Estimation Under Quasi-Experimental Setting');disp(' ');
% Quasi-Experiment
Data_dyadic_PT      = [idx_d, Di Dj A0 A1];
Data_individual_PT  = [idx_i, D, Y0 Y1];

[est_zeta1, est_xi, est_beta_PT, est_pi_PT, est_beta_L_PT, est_beta_S] ...
    = est_PT(Data_dyadic_PT,Data_individual_PT,N);
%est_zeta1
est_xi
est_beta_PT
est_pi_PT
