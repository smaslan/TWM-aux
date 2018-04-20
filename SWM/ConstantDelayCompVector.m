function [CompVector] = ConstantDelayCompVector(fs,FFT_size,delay_s)
% fs [Hz]     : Sampling frequency
% FFT_size    : Size of FFT
% delay_s[sec]:  Cabel-delay to be compensated: in seconds 
% Pos. value remove delay equal to argument, negative will add a delay 

Phi=[];
FFT_half=FFT_size/2;    
    
for x = 0:FFT_size-1
    if x>=FFT_half, Xi=x-FFT_size; else Xi = x; end
    freq=Xi/FFT_size*fs;
    phi=2*pi*freq*delay_s;
    Phi = [Phi phi];
end

CompVector = complex(cos(Phi),sin(Phi));  

