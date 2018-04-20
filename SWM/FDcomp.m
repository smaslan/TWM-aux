
function [ Uo ] = FDcomp(Ui,PhicompVector) % Compensation of Frequenzy-dependent gain and phase errors
% Raw time-frequenzy-time domain. Windowing is done Pre. and Post. this function
% ------------------------------------------
% Ui1: Sampled data for U1 (Voltage channel)
% Ui2: Sampled data for U2 (Current channel)
% PhicompVecto: Complex array, containing phase and gain correction over the spechter

% Check Input arguments
if ~isequal(size(Ui,1), size(PhicompVector,1) , 1)  
    error("In function 'FDcomp(Ui,PhicompVector)', Vector(s) not 1-domentional!");
end
if ~isequal(size(Ui,2), size(PhicompVector,2))  
    error("In function 'FDcomp(Ui1,PhicompVector)', function was called with vector length mismatch!");
end

fft_size=size(Ui,2);

% Check Input array-length
if ~isequal(fft_size, 2^nextpow2(fft_size))  
    error("In function 'FDcomp(Ui1,Ui2,PhicompVector)', Array lengt not exactly 2^N. Length mismatch!");
end

F_domain = fft(Ui,fft_size);  % Time-to-Frequenzy of sampled chunk Current-channel


% Angel = mod(angle(F_domain)+2*pi(),pi());

Fcmpt = F_domain.* PhicompVector;  % Phase and Gain is corrected according to PhiCompVector

Fcmpt(fft_size/2+1) = real(Fcmpt(fft_size/2+1));

%###note: changed for Octave compatibility
%Uo = ifft(Fcmpt,'symmetric'); % Uo2=ifft(YF,fft_size);
Uo = real(ifft(Fcmpt)); % Uo2=ifft(YF,fft_size);


