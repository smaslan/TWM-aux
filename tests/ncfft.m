clc;
clear all;

% test of FFT for non-coherent sampling
% testing scalloping effect

fs = 10000;
N = 1e4;

Ak = 0.5;

fn = linspace(40,60,100);

na = [];
ru = [];
rp = [];
for k = 1:numel(fn)

    f0 = fn(k);
    
    t = [];
    t(:,1) = [1:N-1]/fs;
    
    u1 = sin(t*2*pi*f0);    
    u2 = sin(t*2*pi*f0)*Ak;
    
    [fx,U1,p1] = ampphspectrum(u1,fs,0,0,'flattop_matlab',[],0);
    [fx,U2,p2] = ampphspectrum(u2,fs,0,0,'flattop_matlab',[],0);
    
    [v,id] = min(abs(fx - f0));     
    
    na(end+1) = U2(id);
    ru(end+1) = U2(id)/U1(id);
    rp(end+1) = p2(id) - p1(id);

end

figure
plot(na)
figure
plot(ru)
figure
plot(rp)






