clc;
clear all;

% some testing of modulated signal generation...

fs = 10000;
N = 5000;

f0 = 50;
fm = 5;

A0 = 1;
Am = 0.2;

phm = 10/180*pi;

u  = mod_synth(fs,N,0, f0,A0,0, fm,Am,phm, 'sine', [1 1.5 0.5],[0.1 0.2 -0.3]*pi);

plot(u);


 



