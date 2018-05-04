clc;
clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

%load('wrms_single_tone_bits.matsc','res','vr','p','s')
%load('wrms_spurr_test.matsc','res','vr','p','s')
%load('fftf_test.matsc','res','vr','p','s')
load('fftf_test_2.matsc','res','vr','p','s')



