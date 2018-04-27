clc;
clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

% -----------------------------------------------------------------------------------------------------------------------
% Sensitivity analysis for JV's FFT filter 'sFreqDep_PG_Comp.m'
% -----------------------------------------------------------------------------------------------------------------------

% samples count (2^N):
%  note: 10 means it will simulate random sample counts from (2^9)+1 to 2^10 
s.N_pow = 12;[10:2:20];
% harmonic frequency range relative to fs:
s.f0_rat_min_bin = 6;
s.f0_rat_max = 0.1;[0.02 0.05 0.1 0.2 0.3 0.4 0.45];
% harmonic amplitude:
s.f0_amp = 1;
% randomize phase range [+- rad]:
s.f0_phi_rnd = pi;
% filter max phase [+- rad]:
s.ff_max_phi = 0.1;[0 logspace(-6,-0,10)]*pi;
% filter max real amp dev [+- V/V]:
s.ff_max_amp = 1e-6;[0 logspace(-6,log10(0.5),10)]*pi;
% bits count per f0_amp:
s.bits = 24;[16:2:24];
% rms noise:
s.rms_noise = 1e-6;
% tests per setting:
s.cycles = 1;

proc_FFTF(s)



