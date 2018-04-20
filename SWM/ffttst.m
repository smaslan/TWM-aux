clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

fs = 10;
N = 10;

f = [5]*fs/10;
phi = 0.25*pi;
A = 1;
dc = 0;


t(:,1) = [0:N-1]/fs;

u = dc + sum(A.*sin(t*2*pi*f + phi),2);


U = fft(u)/N
    
Ui = 2*j*fft(u)/N

y1 = U(N/2+1);
y2 = Ui(N/2+1);






return




f0(1,:) = logspace(log10(50.3),log10(2000.3),M);
A = 1e-0*ones(size(f0));
A(1) = 1;
phi = 0;


fm = 0.8;

f_env = 0.5*[0:N-1]/N*fs;
g_env = cos(f_env/fs*2*pi*fm);
p_env = linspace(0.1*pi,0.1*pi,N);

f_gain = interp1(f_env,g_env,f0,'pchip','extrap');
f_phi  = interp1(f_env,p_env,f0,'pchip','extrap');

t(:,1) = [0:N-1]/fs;

u = A.*sin(t.*2.*pi.*f0);
u = sum(u,2);

um = f_gain.*A.*sin(t.*2.*pi.*f0 + f_phi);
um = sum(um,2);


fft_size = 4096;

fft_half = fft_size/2;

fr = [0:fft_half]/fft_size;
ff = [];
fg = interp1(f_env,g_env,fr*fs,'pchip','extrap');
fp = interp1(f_env,p_env,fr*fs,'pchip','extrap');
ff(:,1) = fg.*exp(j*fp);
ff(fft_half+2:fft_size) = conj(ff(end-1:-1:2));
ff(fft_half+1) = ff(fft_half+1);

fc = ff.';

%fc = NullCompVector(fft_size);
%fc = H3458ACompVector(fs,fft_size,100e-6);

[uf,uf,first,last] = sFreqDep_PG_Comp(u',u',fft_size,fc,fc);
uf = uf';

um = um(first:last);

%plot(um(first:last))
%hold on;
%plot(uf,'r')
%hold off;

plot(uf - um)



