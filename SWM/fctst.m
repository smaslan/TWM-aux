
mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

fs = 10000;
N = 20000;

f0 = linspace(5.3,2000.3,N/2);
f0(N/2+1:N) = f0(end:-1:1);
%f0(N/2+1:N) = f0(1:end);
A = 1.0;
phi = 0;

p.fs = fs;
p.jitter = 0;
p.N = N+1;
p.n = 1;
p.Fs = 50;


% - Signal parameters
x = struct(...
 'f0',     5.3,...                                                   % Nominal power frequency
 'harm',   struct('rank',[1],'mag',[A],'pha',[phi]),...                  % Harmonics, phase offset is because one sample is removed later 
 'inter',  struct('freq',[],'mag',[],'pha',[]),...                  % Interharmonics
 'noise',  0,...                                                    % Noise in signal
 %'fMod',   struct('type','triangle','df',df,'f',1/(2*Tramp),'dt',10,'ncycle',2),...     % Frequency modulation
 'fMod',   struct('type','triangle',... % triangle frequency modulation
        'df',2000.3 - 5.3,...                % change of frequency during whole modulation
        'f',0.5/(N/fs),...           % frequency of modulation - function create triangle with -^v-, therefore divide by two, and Tramp is time of half of triangle, therefore again divide by two
        'dt',0,...                      % no effect?
        'ncycle',1),...                 % number of modulation cycles
 'phaMod', struct('type','','dpha',0,'f',0,'dt',0,'ncycle',0),...   % Phase modulation
 'magMod', struct('type','','dmag',0,'f',0,'dt',0,'ncycle',0));     % Amplitude modulation
%

u= Waveform_Generator2(p,x);
u = u.y(2:end);

     
t = [0:N-1]/fs;

%u = A*sin(t.*2.*pi.*f0 + phi);

Am = 0.5;
fm = 20;

um = u.*(1 + Am*sin(f0/(fs)*2*pi*fm));


fft_size = 2048;

fft_half = fft_size/2;

fr = [0:fft_half]/fft_size;

ff = [];
ff(:,1) = sin(fr*2*pi*fm)*Am + A;
ff(fft_half+2:fft_size) = ff(end-1:-1:2);

fc = ff.';



%fc = NullCompVector(fft_size);
%fc = H3458ACompVector(fs,fft_size,100e-6);

[uf,uf,first,last] = sFreqDep_PG_Comp(u,u,fft_size,fc,fc);

plot(um(first:last))
hold on;
plot(uf,'r')
hold off;



