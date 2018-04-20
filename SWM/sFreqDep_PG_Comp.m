% Freqency Dependant Phase and Gain Compensation 
function [Uc1, Uc2, first, last] = sFreqDep_PG_Comp(U1,U2,fft_size,CmpVector1,CmpVector2)
% Input arguments:
% U1,U2        : Sampling buffer data arrays (uncompensated) 
% fft_size     : FFT size.
% CmpVector1,CmpVector2: Complex vector holding Gain and Phase compensation

% Output: 
% Uc1, Uc2    : (Compensated output data arrays.)
% first, last : Index of which part of the input buffer is used for output

if ~isequal(fft_size, 2^nextpow2(fft_size))  
    error("In function 'FDcomp(Ui1,Ui2,PhicompVector)', Array lengt not exactly 2^N. Length mismatch!");
end
N = 2^nextpow2(fft_size); % Size of FFT
Hw = HannFcompMask(N)/2;    % Masking for overlapping processing windows



% set startconditions
C1=zeros(1,N/4); 
C2=zeros(1,N/4); 
Uc1=[]; % Start with empty output array
Uc2=[]; % Start with empty output array
[Bs Be Tp]=PackMan(size(U1,2),N,N/4,1);           
%[Bs Be Tp]
for frame=1:Tp
    % find frame:
    [Bs Be Tp]=PackMan(size(U1,2),N,N/4,frame);  %Find Frame Position         
    %[frame Bs Be]
    
    Uo1 = FDcomp(U1(Bs:Be),CmpVector1); % Run compensation over frame data:
    Uo2 = FDcomp(U2(Bs:Be),CmpVector2); % Run compensation over frame data:
    % Windowing the result and 
    Uo1= Uo1.*Hw; Uo2= Uo2.*Hw;  
    
    if frame ~=1,
        % pick data for sum and concatenation
        Uc1=[Uc1,C1+Uo1(N*1/4+1:N*2/4)]; %data for morhping with prew. frame
        Uc2=[Uc2,C2+Uo2(N*1/4+1:N*2/4)]; %data for morhping with prew. frame
    end
    C1=Uo1(N*2/4+1:N*3/4); %data for morhping with next frame
    C2=Uo2(N*2/4+1:N*3/4); %data for morhping with next frame
end

first =N/2+1;
last = Be-N/2;


