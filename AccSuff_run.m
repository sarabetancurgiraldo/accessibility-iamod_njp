close all; clear; clc;

pwd
addpath(genpath("/home/mech001/20222295/Sara/YALMIP-master"));
fclose('all');

load('model/data_g.mat');
load("model/data_shortPaths.mat");

%% Variables

nCarRange = [0 2e3 3e3 4e3 5e3];
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

AccSuff = cell(Ns);%(Ts,nC,3);

for i_nCar = 1:nC
nCar = nCarRange(i_nCar); 

for i_Tsuff = 1:Ts
Tsuff = TsuffRange(i_Tsuff);

%% Optimizations

if nCar ~= 0

for i_Nsuff = 1:Ns
Nsuff = NsuffRange(i_Nsuff);

% AccessibilitySufficiency
str_save = sprintf('output/Nsuff/%d/nCar/%d/Tsuff/%d/AccSuff.mat',Nsuff,nCar,Tsuff*60);
AccessibilitySufficiency(G,B,pc_unique,X_Iamod,X_amod,D,nOD,R_selector,Nsuff, ...
                         population_region,str_save,alpha,Tsuff,nArcs,nCar)

%% Path flow allocation

% AccessibilitySufficiency
load(sprintf('output/Nsuff/%d/nCar/%d/Tsuff/%d/AccSuff.mat',Nsuff,nCar,Tsuff*60));
X_matrix = sol_AccSuff.X;
t_AvgAcc = (sol_AccSuff.F_Iamod.*sol_AccSuff.t_Iamod + ...
            sol_AccSuff.F_amod.*sol_AccSuff.t_amod)./...
            (sol_AccSuff.F_Iamod+sol_AccSuff.F_amod);
epsilonAcc = (max(0,t_AvgAcc-Tsuff)/Tsuff).^2;
epsilonAcc = epsilonAcc(1:nOD);
fp_save = sprintf('output/Nsuff/%d/nCar/%d/Tsuff/%d/path_flows_AccSuff.mat',Nsuff,nCar,Tsuff*60);
% [AccSuff(1,1,1),...
%  AccSuff(1,1,2),...
%  AccSuff(1,1,3)] = path_flows_Leo(Tsuff,X_matrix,epsilonAcc',fp_save,D,B,G);
[AccSuff{i_Nsuff,i_Tsuff,i_nCar,1},...
 AccSuff{i_Nsuff,i_Tsuff,i_nCar,2},...
 AccSuff{i_Nsuff,i_Tsuff,i_nCar,3}] = path_flows_Leo(Tsuff,X_matrix,epsilonAcc',fp_save,D,B,G);

% str_save = sprintf('output/Nsuff/%d/nCar/%d/Tsuff/%d/J.mat',Nsuff,nCar,Tsuff*60);
% save(str_save,'AccSuff');

%% Heatmap

% AccessibilitySufficiency
load(sprintf('output/Nsuff/%d/nCar/%d/Tsuff/%d/AccSuff.mat',Nsuff,nCar,Tsuff*60))
E_Iamod = max(0,60*(sol_AccSuff.t_Iamod-Tsuff)).^2;
E_amod = max(0,60*(sol_AccSuff.t_amod-Tsuff)).^2;
epsilonAcc = (sol_AccSuff.F_Iamod.*E_Iamod + ...
              sol_AccSuff.F_amod.*E_amod)/sum(alpha);
fp_load = sprintf('output/Nsuff/%d/nCar/%d/Tsuff/%d/path_flows_AccSuff.mat',Nsuff,nCar,Tsuff*60);
fp_save = sprintf('output/Nsuff/%d/nCar/%d/Tsuff/%d/AFI_heatmap_AccSuff.mat',Nsuff,nCar,Tsuff*60);
AFI_heatmap_sq(Tsuff,fp_load,fp_save,epsilonAcc',D,true)

end
end
end
end

str_save = sprintf('output/J_AccSuff.mat');
save(str_save,'AccSuff');
