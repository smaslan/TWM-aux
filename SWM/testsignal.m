
function s=testsignal(N,fs,bf,amp,dc,phi,noiseamp)
% N,fs,bf,amp,phi,noiseamp
%
% N        : Numbers of samples
% fs       : Sampling frequency
% bf       : Base signal frequency
% amp      : Base signal Amplitude (peek)
% dc       : Base signal DC-offset
% phi      : Base signal phase in degrees
% noiseamp : random noise amplitude
 a=amp*sqrt(2);
 k=2*pi*bf/fs;
 p=2*pi/360*phi;
for i = 1:N
    s(i)=a*sin(k*i+p)+(rand()-0.5)*noiseamp+dc;
end
