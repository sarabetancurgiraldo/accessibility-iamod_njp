close all; clear; clc;

pwd
addpath(genpath("/home/mech001/20222295/Sara/YALMIP-master"));
fclose('all');

load('model/data_g.mat');
load("model/data_shortPaths.mat");

%% Variables

nCarRange = [0 3e3 4e3 5e3];
TsuffRange = [15/60 20/60 25/60];
NsuffRange = [30 35 40];

% alpha matrix - # of trips per hour for each od-pair 
alpha           = sum(abs(D),1)/2;
t               = G.Edges.Weight;
X_Iamod         = Xfast;
X_amod          = Xslow;

% nOD             = 5;
% D               = D(:,1:nOD);
% X_Iamod         = X_Iamod(:,1:nOD);
% X_amod          = X_amod(:,1:nOD);
% R_selector      = R_selector(:,1:nOD);
% alpha           = alpha(:,1:nOD);

%% Create object with data
nC = length(nCarRange);
Ts = length(TsuffRange);
Ns = length(NsuffRange);

UtilEff = zeros(Ts,nC,3);
CommSuff = zeros(Ts,nC,3);
TripSuff = zeros(Ts,nC,3);
AccSuff = zeros(Ts,nC,3);

nCar = 4000;
Tsuff = 15/60;
i_Tsuff = 1;
i_nCar = 3;

% CommuteSufficiency
load(sprintf('output/nCar/%d/Tsuff/%d/CommSuff.mat',nCar,Tsuff*60));
X_matrix = sol_comSuff.X;
X_matrix = X_matrix(:,1:nOD);
epsilonComm = sol_comSuff.epsilon;
epsilonComm = epsilonComm(1:nOD);
fp_save = sprintf('output/nCar/%d/Tsuff/%d/path_flows_CommSuff.mat',nCar,Tsuff*60);
[CommSuff(i_Tsuff,i_nCar,1), ...
 CommSuff(i_Tsuff,i_nCar,2), ...
 CommSuff(i_Tsuff,i_nCar,3)] = path_flows_Leo(Tsuff,X_matrix,epsilonComm,fp_save,D,B,G);


