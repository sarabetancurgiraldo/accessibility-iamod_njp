close all; clear; clc;

pwd
addpath(genpath("/home/mech001/20222295/Sara/YALMIP-master"));
fclose('all');

load('model/data_g.mat');
load("model/data_shortPaths.mat");

%% Variables

nCarRange = [4e3];
TsuffRange = [20/60];
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

for i_nCar = 1:nC
nCar = nCarRange(i_nCar); 

for i_Tsuff = 1:Ts
Tsuff = TsuffRange(i_Tsuff);

%% Optimizations

% UtilitarianEfficiency
str_save = sprintf('output/nCar/%d/Tsuff/%d/UtilEff.mat',nCar,Tsuff*60);
UtilitarianEfficiency(nCar,G,B,nArcs,D,str_save);

%% Path flow allocation

% UtilitarianEfficiency
load(sprintf('output/nCar/%d/Tsuff/%d/UtilEff.mat',nCar,Tsuff*60));
X_matrix        = sol_utilEff.X;
X_matrix        = X_matrix(:,1:nOD);
epsilonUEff     = (max(0,(2*(t'*X_matrix)'./sum(abs(D),1)')-Tsuff)/Tsuff).^2;
fp_save = sprintf('output/nCar/%d/Tsuff/%d/path_flows_UtilEff.mat',nCar,Tsuff*60);
[UtilEff(i_Tsuff,i_nCar,1), ...
 UtilEff(i_Tsuff,i_nCar,2), ...
 UtilEff(i_Tsuff,i_nCar,3)] = path_flows_Leo(Tsuff,X_matrix,epsilonUEff,fp_save,D,B,G);

%% Heatmap

% UtilitarianEfficiency
load(sprintf('output/nCar/%d/Tsuff/%d/UtilEff.mat',nCar,Tsuff*60));
X_matrix = sol_utilEff.X;
tm_UtilEff = (2*(t'*X_matrix)'./sum(abs(D),1)');
epsilonUEff = (max(0,60*(tm_UtilEff-Tsuff))).^2;
fp_load = sprintf('output/nCar/%d/Tsuff/%d/path_flows_UtilEff.mat',nCar,Tsuff*60);
fp_save = sprintf('output/nCar/%d/Tsuff/%d/AFI_heatmap_UtilEff.mat',nCar,Tsuff*60);
AFI_heatmap_sq(Tsuff,fp_load,fp_save,epsilonUEff,D,false)

end
end

str_save = sprintf('output/nCar/%d/J_UE.mat',nCar);
save(str_save,'UtilEff');%,'CommSuff','TripSuff');%,'AccSuff');


