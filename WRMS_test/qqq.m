clc;
clear all;
close all;


fs = 100000;

N = 100000;

t(:,1) = [0:N-1]/N*fs;

tic;
for k = 1:100
    f  = rand(1,10);
    A  = rand(1,10);
    ph = rand(1,10)*2*pi;
    
    w = 2*pi*f;
    u = sum(A.*sin(w.*t + ph),2);
end
toc()


return





rms_noise = 1e-6;

N = 100000;

tw(:,1) = [0:N-1]/N*2*pi;

u = rms_noise*randn(N,1);
%u = 2^0.5*sin(10*tw);

w = window_coeff('blackman',N);

% noise rms in time domain:
m_noise_t = mean((u.*w).^2)^0.5/mean(w.^2)^0.5




% some badass window for FFT:
w = window_coeff('flattop_248D',N)(:);
wg = mean(w);


% FFT:
U = fft(u.*w)/N*2;
U = abs(U(1:floor(N/2)));
U = U/wg;

% RMS noise from windowed spectrum:
m_noise_f = sum(0.5*(U).^2)^0.5/mean(w.^2)^0.5*wg



