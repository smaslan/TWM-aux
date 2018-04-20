function  maskWindow=HannFcompMask(N) % N = fft data buffer size, W is Width of Han-win. (Default N/4)

NFFTh = 2^nextpow2(N)/4; % Next power of 2 from length of y 
maskWindow= [zeros(1,NFFTh) hanningw(NFFTh*2) zeros(1,NFFTh)];
 
