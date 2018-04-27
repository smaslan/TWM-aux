clc;
clear all;
close all;



rms_noise = 1e-6;

N = 100000;

tw(:,1) = [0:N-1]/N*2*pi;

u = rms_noise*randn(N,1);
i = rms_noise*randn(N,1);
%u = 2^0.5*sin(10*tw);

w = window_coeff('blackman',N);
W = mean(w.^2)^-0.5;

% noise rms in time domain:
m_noise_t = mean(u.*i.*w.^2)*W^2


U = fft(u)/N*2;
I = fft(i)/N*2;

sum(0.5*U.*I).^0.5


return

% some badass window for FFT:
w = window_coeff('flattop_248D',N)(:);
wg = mean(w);


% FFT:
U = fft(u.*w)/N*2;
U = abs(U(1:floor(N/2)));
U = U/wg;

% RMS noise from windowed spectrum:
m_noise_f = sum(0.5*(U).^2)^0.5/mean(w.^2)^0.5*wg



