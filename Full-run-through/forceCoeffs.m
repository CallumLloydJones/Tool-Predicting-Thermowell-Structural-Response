clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data

mesh= ["Straight\Re1e5" "Helical\Re1e5"];
U= [3.5 3.5];
D = 30e-3;
start = [3.0 3.0];

St= [];
Re = [1e3 1e4 1e5];

dir = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\';
rootdir = dir+mesh;
for i = 1:length(mesh)
    rootdirs = dir+mesh(i);
    rest = plotliftdrag(rootdirs, mesh(i), start(i), U(i), D);
    %St = [St rest];
end


