clear all

baseroot = pwd();
chdir(baseroot);
addpath('mfiles');

Dyadic_data_raw = importdata([baseroot,'\Dyadic_data.xlsx']);
Ind_data_raw = importdata([baseroot,'\Ind_data.xlsx']);

Data_dyadic     = Dyadic_data_raw.data;
Data_individual = Ind_data_raw.data;
dyadic_index    = Data_dyadic(:,2)~=Data_dyadic(:,3);
Data_dyadic     = Data_dyadic(dyadic_index,:);

% Load from dyadic data
idx_d   = Data_dyadic(:,1:3);
A0      = Data_dyadic(:,4);
A1      = Data_dyadic(:,5);
DA      = A1-A0;
Di      = Data_dyadic(:,6);
Dj      = Data_dyadic(:,7);
Gdum_d  = Data_dyadic(:,8:26);

% Load from individual data
idx_i = Data_individual(:,[1,2]);
Y0      = Data_individual(:,3);
Y1      = Data_individual(:,4);
DY      = Y1-Y0;
D       = Data_individual(:,5);
Gdum_i  = Data_individual(:,6:23);

%%
clc
N = 20;
Data_dyadic_RE      = [idx_d, Di, Dj, A1];
Data_individual_RE  = [idx_i, D, Y1];
[est_zeta1, est_beta_RE, est_pi_RE, est_beta_L_RE, est_beta_S] ...
    = est_RE(Data_dyadic_RE,Data_individual_RE,[Gdum_i Y1>0],N,"g");
est_zeta1
est_beta_RE
est_pi_RE


Data_individual_RE  = [idx_i, D, log(Y1+1)];
[est_zeta1, est_beta_RE, est_pi_RE, est_beta_L_RE, est_beta_S] ...
    = est_RE(Data_dyadic_RE,Data_individual_RE,[Gdum_i Y1>0],N,"g");
est_beta_RE
est_pi_RE