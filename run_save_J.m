close all; clear; clc;

% pwd
% addpath(genpath("/home/mech001/20222295/Sara/YALMIP-master"));
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


for i_nCar = 1:nC
nCar = nCarRange(i_nCar); 

for i_Tsuff = 1:Ts
Tsuff = TsuffRange(i_Tsuff);

%% Path flow allocation

% UtilitarianEfficiency
load(sprintf('output/nCar/%d/Tsuff/%d/UtilEff.mat',nCar,Tsuff*60));
X_matrix        = sol_utilEff.X;
X_matrix        = X_matrix(:,1:nOD);
epsilonUEff     = (max(0,(2*(t'*X_matrix)'./sum(abs(D),1)')-Tsuff)/Tsuff).^2;
fp_save = sprintf('output/nCar/%d/Tsuff/%d/path_flows_UtilEff.mat',nCar,Tsuff*60);
[UtilEff(i_Tsuff,i_nCar,1), ...
 UtilEff(i_Tsuff,i_nCar,2), ...
 UtilEff(i_Tsuff,i_nCar,3)] = path_flows_values(Tsuff,X_matrix,epsilonUEff,fp_save,t,D);

if nCar ~= 0

% CommuteSufficiency
load(sprintf('output/nCar/%d/Tsuff/%d/CommSuff.mat',nCar,Tsuff*60));
X_matrix = sol_comSuff.X;
X_matrix = X_matrix(:,1:nOD);
epsilonComm = sol_comSuff.epsilon;
epsilonComm = epsilonComm(1:nOD);
fp_save = sprintf('output/nCar/%d/Tsuff/%d/path_flows_CommSuff.mat',nCar,Tsuff*60);
[CommSuff(i_Tsuff,i_nCar,1), ...
 CommSuff(i_Tsuff,i_nCar,2), ...
 CommSuff(i_Tsuff,i_nCar,3)] = path_flows_values(Tsuff,X_matrix,epsilonComm,fp_save,t,D);

% TripSufficiency 
load(sprintf('output/nCar/%d/Tsuff/%d/TripSuff.mat',nCar,Tsuff*60))
X_matrix = sol_Tripsuff.X;
X_matrix = X_matrix(:,1:nOD);
t_AvgTrip = (sol_Tripsuff.F_Iamod.*sol_Tripsuff.t_Iamod + ...
             sol_Tripsuff.F_amod.*sol_Tripsuff.t_amod)./...
             (sol_Tripsuff.F_Iamod+sol_Tripsuff.F_amod);
epsilonTrip = (max(0,t_AvgTrip-Tsuff)/Tsuff).^2;
epsilonTrip = epsilonTrip(1:nOD);
fp_save = sprintf('output/nCar/%d/Tsuff/%d/path_flows_TripSuff.mat',nCar,Tsuff*60);
[TripSuff(i_Tsuff,i_nCar,1), ...
 TripSuff(i_Tsuff,i_nCar,2), ...
 TripSuff(i_Tsuff,i_nCar,3)] = path_flows_values(Tsuff,X_matrix,epsilonTrip',fp_save,t,D);

end
end
end


str_save = sprintf('output/J.mat');
save(str_save,'UtilEff','CommSuff','TripSuff');




function [TT,AFI,AFImean] = path_flows_values(T_max,X_matrix,epsilon,fp_save,t,D)

    load(fp_save,'t_paths_cell','f_paths_cell')

    minAcc = zeros(size(X_matrix,2),1);
    
    for od_pair_idx = 1:size(X_matrix,2) 
    
        t_paths = t_paths_cell{od_pair_idx,1};
        minAcc_f = f_paths_cell{od_pair_idx,1};
        
        dt_max0 = t_paths-T_max;
        dt_max0(dt_max0<0) = 0;
        minAcc(od_pair_idx) = (dt_max0/T_max).^2'*minAcc_f;

    end

    TT = 60*t'*X_matrix*ones(size(D,2),1)/(sum(sum(abs(D)))/2);
    AFI = 60*sum(minAcc)/(sum(abs(D),'all')/2);
    AFImean = 60*epsilon'*sum(abs(D),1)'/(sum(abs(D),'all'));

end