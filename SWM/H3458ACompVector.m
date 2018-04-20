function [CompVector] = H3458ACompVector(fs,FFT_size,intgration_time)
%Calculate Compensation for 3458A-Frequency-dependent Gain:
% fs [Hz] : Sampling frequency
% FFT_size: Size of FFT
% intgration_time [sec]:  HP3458A Apperture time: example: intgration_time=90us 
    
Cinv=[];
FFT_half=FFT_size/2;    
    
for x = 0:FFT_size-1
    if x>=FFT_half, Xi=FFT_size-x; else Xi = x; end
    q = Xi/FFT_size*pi()*fs*intgration_time;
    cinv=1/sincm(q);
    Cinv=[Cinv cinv];
end
plot(Cinv)
CompVector=complex(Cinv,zeros(1,FFT_size));


