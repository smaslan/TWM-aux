clc;
clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
%addpath([mfld filesep() 'var']);



rms_noise = 1e-3;

ab_rat = logspace(-8,0,100);

A = 1;
B = A*ab_rat;

r = ((rms_noise/A)^2 + (rms_noise./B).^2).^0.5;

loglog(ab_rat,r)



