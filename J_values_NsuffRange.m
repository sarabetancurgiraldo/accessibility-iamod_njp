close all; clear; clc;

load('model/data_g.mat');
load("model/data_shortPaths.mat");

nCar = 4000;
Tsuff = 20/60;
NsuffRange = [30 35 40];
Ns = length(NsuffRange);

for i_Nsuff = 1:Ns
Nsuff = NsuffRange(i_Nsuff);


load(sprintf('output/Nsuff/%d/nCar/%d/Tsuff/%d/J.mat',Nsuff,nCar,Tsuff*60));
fp_load = sprintf('output/Nsuff/%d/nCar/%d/Tsuff/%d/AccSuff.mat',Nsuff,nCar,Tsuff*60);
load(fp_load);
load(sprintf('output/Nsuff/%d/nCar/%d/Tsuff/%d/AFI_heatmap_AccSuff.mat',Nsuff,nCar,Tsuff*60));

Tavg = AccSuff(1,1,1); 

sprintf(['AccSuff: \n' ...
        'Travel time: %.4f'], Tavg)

pause
end