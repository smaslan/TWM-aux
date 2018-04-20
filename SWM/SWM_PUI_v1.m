% SWM Example code
% 
% 1) Data_in_func(); % Two arrays of sampled or simulated 3458A voltage-data
%    Simulated data used in this example: testsignal(N,fs,Hz,U,Udc,dPhi,noise);

% 2) Data has frequency-dependent errors (phase-Gain), compensation by function 
%     sFreqDep_PG_Comp(U1,U2,fft_size,CmpVector1,CmpVector2);
%     NB: Can only compensate data between Data(fft_size/2:SIZE-~fft_size/2) 
%     The function returns two compensated arrays and the index of START and END.
%     In this example, compensation-arrays simulates cabel delay: 
%     func: linearPhaseCompVector(fft_size,rad);

% 3) Timebased calculation of Base Properties; 
%    [U1rmws U2rmws U1dc U2dc Pact Prea Papp PF] = swm1(U1,U2,rootedWW,WW)
%    
%----Outputs: -----------------
%    U1rmws:  RMS-value of U1(Voltage-channel)
%    U2rmws:  RMS-value of U1(Voltage-channel)
%    U1dc:    DC-offset of U1 (Voltage-channel)
%    U2dc:    DC-offset of U2 (Current-channel)
%    Pact:    RMS-value of Active Power
%    Prea:    RMS-value of Reactive Power
%    Papp:    RMS-value of Apparent Power
%    PF:      Power factor
%
%  ACCURRACY: for better then max. 1ppm error, the input data length must
%  be longer then 31.75 periodes of the input signal Base Frequency, for
%
%  example:
%  50Hz, sampled at 10kHz => min. 6350 samples, in addition add ca.
%  1*fft-size since the compensates signal is shorter then the original
%  signal


% Simulated Signal Properties 
fs = 10000;%10000;     % Sampling frequenzy

N=(fs/50*32.0+2048);  % Sample-length(number of sampels) For 50Hz 10kHz clock, N more then 6400 > ca. 1ppm max calc-error.
Hz=50.010273;      % Hz (frequenzy of simulated input signal)
U=230.0*1;        % Volt RMS of simulated input signal
I=230.00*1;        % Amp. RMS of simulated input signal
Udc= 0.0;       % DC-Volt of simulated input signal
Idc= -0.0;       % DC-Amp. of simulated input signal
phi=60.0;       % Phi. angel[deg] between Voltage and Current. PF= Ind. 0.87 of simulated input signal

% seting up simulated signals
dPhi= rand()*360;    % for simulation: add randomness to start-phase
U1=testsignal(N,fs,Hz,U,Udc,dPhi,0);      % Generating U1 (simulated Signal-array for U1 (Voltage-channel)
U2=testsignal(N,fs,Hz,I,Idc,dPhi+phi,0);  % Generating U1 (simulated Signal array for U2 (Current-channel)

%  Compensation of Freqency Dependant Phase & Gain (U1(1:2048),U2(1:2048),Cphi,CGain)
fft_size=1024;
fft_size=2^nextpow2(fft_size);  %Ensure 2^N - length
fft_overlapCrop = fft_size/4;

% Examples of Frequenzy dependent compensation of Phase and Gain;

% Null compensation: No change in Phase and Gain;  Simulated input: 
CmpVector1 = NullCompVector(fft_size); % No compensation for Channel 1
CmpVector2 = NullCompVector(fft_size); % No compensation for Channel 2

% Cable Delay Compensation:  ConstantDelayCompVector(fs,fft_size,delay[sek]); 
%CmpVector1 = ConstantDelayCompVector(fs,fft_size,10e-7); % Compensation for Cable delay Channel 1
%CmpVector2 = ConstantDelayCompVector(fs,fft_size,10e-7); % Compensation for Cable delay Channel 2

% Compensation for 3458A-Frequency-dependent Gain: H3458ACompVector(fs,FFT_size,intgration_time)
intgration_time=90e-6; % HP3458A Apperture time: example intgration_time=90us 
%CmpVector1=H3458ACompVector(fs,fft_size,90e-6); % Gain and Phase compensation vector Channel 1
%CmpVector2=H3458ACompVector(fs,fft_size,90e-6); % Gain and Phase compensation vector Channel 1


%Call Compensation function with correct compensation vector;
[Uc1, Uc2, first, last] = sFreqDep_PG_Comp(U1,U2,fft_size,CmpVector1,CmpVector2);
                                         %Uc1=U1; Uc2=U2;last=size(U1,2);first=1; 

% one time window generation (can be user repetedly, not calculate it every sycle) 
WW= hanningw(last-first+1);  rootedWW = sqrt(WW);  %SumWW=sum(WW);

% ===== SWM: Calculation of power Parameters ==================
[U1rmws U2rmws U1dc U2dc Pact Prea Papp PF] = swm1(Uc1, Uc2, rootedWW, WW);
% ===== SWM ===================================================

% Test of Result compared with theoretically correct values. 

% Calculate Teoretical correct value: (for comparison with calculated results)
cUrms=sqrt(Udc^2+(U*sqrt(2))^2/2); % Comb. Ac and Dc Volt RMS  
cIrms=sqrt(Idc^2+(I*sqrt(2))^2/2); % Comb. Ac and Dc Amp. RMS  
%[cUrms cIrms]
cPact= U*I*cos(phi/180*pi)+Udc*Idc;% Active Power
cPrea= U*I*sin(phi/180*pi);                        % Reactive Power
cPapp= cUrms*cIrms;                                % Apparent Power

eU1 =(U1rmws-cUrms)/cUrms*1000000 ;       % Error U1 calc [ppm]
eU2 =(U2rmws-cIrms)/cIrms*1000000 ;       % Error U2 calc [ppm]

% U=Udc; I=Idc; % Test DC
Papp2= sqrt(Pact^2+Prea^2)
ePact  = (Pact-cPact)/cPapp*1000000;% Error of Active power calculation[ppm]
ePrea  = (Prea-cPrea)/cPapp*1000000;% Error of Reavtive power calculation[ppm]
ePapp  = (Papp2-Papp)/Papp*1000000; % Error of Reactive Power calculation[ppm]

% print result values
U1rmws, U2rmws, U1dc, U2dc, Pact, Prea, Papp, PF, eU1, eU2, ePact, ePrea, ePapp % mean((U1(first:last).*WW'))-Udc, mean((U2(first:last).*WW'))-Idc

