close all;
clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

fs = 10000;
N = 10000;
M = 10;

f0(1,:) = logspace(log10(50.3),log10(2000.3),M);
A = 1e-0*ones(size(f0));
A(1) = 1;
phi = rand(size(A))*2*pi;
dc = 0.001;


fm = 0.0;

f_env = 0.5*[0:N]/N*fs;
g_env = 1*cos(f_env/fs*2*pi*fm);
p_env = linspace(0.0*pi,0.001*pi,N+1);

f_gain = interp1(f_env,g_env,f0,'pchip');
f_phi  = interp1(f_env,p_env,f0,'pchip');

t(:,1) = [0:N-1]/fs;

u = A.*sin(t.*2.*pi.*f0 + phi);
u = dc + sum(u,2);

um = f_gain.*A.*sin(t.*2.*pi.*f0 + phi + f_phi);
um = dc + sum(um,2);


%fft_size = 4096;
fft_size = 2^nextpow2(N/4)

fft_half = fft_size/2;

fr = [0:fft_half]/fft_size;
ff = [];
fg = interp1(f_env,g_env,fr*fs,'pchip');
fp = interp1(f_env,p_env,fr*fs,'pchip');
NF = round(fft_size*0.01)
fp(end-NF:end) = fp(end-NF:end).*(0.5 + 0.5*cos([0:NF]/NF*pi));
ff(:,1) = fg.*exp(j*fp);
ff(fft_half+2:fft_size) = conj(ff(end-1:-1:2));
ff(1) = fg(1);
%ff(fft_half+1) = ff(fft_half+1)*2;
ff(fft_half+1) = 1*fg(end)*cos(fp(end)); % bad approximation
%ff(fft_half+1) = (ff(fft_half+1) + conj(ff(fft_half+1)))*1;

fc = ff.';

%fc = NullCompVector(fft_size);
%fc = H3458ACompVector(fs,fft_size,100e-6);

[uf,uf2,first,last] = sFreqDep_PG_Comp(u',u',fft_size,fc,fc);
uf = uf';

um = um(first:last);

%plot(um(first:last))
%hold on;
%plot(uf,'r')
%hold off;

%figure;
%plot(uf - um)


S = numel(uf);
[fc,Ac,phc] = ampphspectrum(uf, fs, 0, 0, 'flattop_248D', 0, 0);
[fc,At,pht] = ampphspectrum(um, fs, 0, 0, 'flattop_248D', 0, 0);

% harmonic frequencies DFT bins:
fid = round(f0/fs*S) + 1;

dA = At(fid)./Ac(fid) - 1
dA = mod(phc(fid) - pht(fid) + pi,2*pi)-pi
 



figure;
loglog(fc,At)
hold on;
plot(fc,Ac,'r')
hold off;




