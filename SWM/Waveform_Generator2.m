function [signal,waveform,ref]=Waveform_Generator2(param,x)
%
%
%   INPUTS:
%   % --- Acquisition parameters
%   param = struct(...
%    'fs',     0,...     % Sampling frequency (ADC)
%    'jitter', 0,...     % Sampling jitter (nS)
%    'N',      0,...     % Total number of samples acquired
%    'n',      0,...     % Number of samples before t0
%    'Fs',     0,...     % PMU frames per second
%    The hyohithesis is that the sampling frequency is a multiple
%    of the PMU frame per seconds
%
% --- Signal parameters
%   x = struct(...
%    'f0',     50,...                                                   % Nominal power frequency
%    'harm',   struct('rank',[],'mag',[],'pha',[]),...                  % Harmonics 
%    'inter',  struct('freq',[],'mag',[],'pha',[]),...                  % Interharmonics
%    'noise',  0,...                                                    % Noise in signal
%    'fMod',   struct('type','','df',0,'f',0,'dt',0,'ncycle',0),...     % Frequency modulation
%    'phaMod', struct('type','','dpha',0,'f',0,'dt',0,'ncycle',0),...   % Phase modulation
%    'magMod', struct('type','','dmag',0,'f',0,'dt',0,'ncycle',0));     % Amplitude modulation
%
%   OUTPUTS
% --- Sampled waveform
%   signal = struct(...
%    'time',[],...       % Time
%    'index',[],...      % Position index) in time array)corresponding to PMU frames (first sample)
%    'y',[]);            % Samples
%
% --- Sampled waveform parameters
%   waveform = struct(...
%    'inst_mag',[],...     % Instantantaneous magnitude of fundamental
%    'inst_pha',[],...     % Instantantaneous phase of fundamental
%    'inst_freq',[],...    % Instantaneous frequency of fundamental
%    'inst_ROCOF',[]);     % Instantaneous rate of frequecy change
%  
% --- Reference values (Theoretical syncrophasor)
%   ref = struct(...
%    'timeF',[],...       % Reporting time of PMU (FS)
%    'index',[],...      % Position index(in signal)of frame samples
%    'mag',[],...        % Magnitude
%    'pha',[],...        % Phase
%    'freq',[],...       % Frequency
%    'ROCOF',[]);        % Rate of frequency change

%   Version 2.0
%   14 November 2011 
%   the phases are in degrees
%
%   last modification 12 December 2011
%   3th Htpothesis with copmuted NW (only as comment) 
%
%   Derived  from the function GenerateSignal.m
%   Jean-Pierre Brown (March 2011)

%
% Evaluation of the signal amd waveform stuctures
% Some name of the variables have been changed
%

% ---Index and times with the additional jitter

signal.index=1:1:param.N;
jit=param.jitter*1e-9*(2*rand(1,max(size(signal.index)))-1.0);
signal.time = (-param.n:1:param.N-param.n-1)*(1/param.fs)+jit;

% --- Frequency of fundamental
waveform.inst_freq = ones(1,param.N)*x.f0;
if ~strcmp(x.fMod.type,'')
    waveform.inst_freq(param.n+1:param.N) = modulate(x.fMod.type,...
                waveform.inst_freq(param.n+1:param.N),x.fMod.df,...
                x.fMod.f,signal.time(param.n+1:param.N),...
                x.fMod.ncycle);
end

% --- Phase of fundamental
signalpha = ones(1,param.N)*x.harm.pha(1); %corrected sync
if ~strcmp(x.phaMod.type,'')
    signalpha(param.n+1:param.N) = modulate(x.phaMod.type,...
                signalpha(param.n+1:param.N),x.phaMod.dpha,...
                x.phaMod.f,signal.time(param.n+1:param.N),...
                x.phaMod.ncycle);
end

% --- Magnitude of fundamental
waveform.inst_mag = ones(1,param.N)*x.harm.mag(1); %corrected sync
if ~strcmp(x.magMod.type,'')
    waveform.inst_mag(param.n+1:param.N) = modulate(x.magMod.type,...
                waveform.inst_mag(param.n+1:param.N),x.magMod.dmag,...
                x.magMod.f,signal.time(param.n+1:param.N),...
                x.magMod.ncycle);
end

% --- Generation of fundamental and harmonics
signal.y = zeros(1,max(size(signal.time)));
for i = 1:max(size(x.harm.rank))           %corrected sync
    if i == 1
        dw = 2*pi*[0 waveform.inst_freq/param.fs];
        dw = dw(1:max(size(dw))-1);
        aw = 2*pi*[0 diff(waveform.inst_freq)]*(1/param.fs)^2;      
        phase = cumsum(dw+aw);
        signal.y = signal.y + waveform.inst_mag.*cos(phase+signalpha*pi/180);
%Instantaneous phase in degree not wrapped        
% ***** previous    waveform.inst_pha = phase+signalpha*2*pi/360;  *****
        waveform.inst_pha = phase*180/pi+signalpha;  
    
    else
        signal.y = signal.y...
                 + x.harm.mag(i)*cos(x.harm.rank(i)*phase...
                 + x.harm.pha(i)*2*pi/360);
    end
end

% --- Generatation of interharmonics
if max(size(x.inter.freq))>= 1
    for i = 1:max(size(x.async.freq))
        signal.y = signal.y...
                 + x.inter.mag(i)*(cos(2*pi*x.inter.freq(i)*signal.time...
                 + x.inter.pha(i)*2*pi/360));
    end
end

% --- Generation of white noise
signal.y = signal.y+x.noise*randn(1,max(size(signal.time)));

%
% Evaluation of the waveform stucture
%

waveform.inst_ROCOF=([waveform.inst_freq(2)-waveform.inst_freq(1) diff(waveform.inst_freq)]+....
    [diff(waveform.inst_freq) waveform.inst_freq(param.N)-waveform.inst_freq(param.N-1)])*param.fs/2 ;
%
% Evaluation of the structure ref
%

signal.index=param.n:(param.fs/param.Fs):param.N;
ref.index=signal.index;
ref.timeF=signal.time(ref.index);
ref.mag=waveform.inst_mag(ref.index);
ref.pha=waveform.inst_pha(ref.index);

% The angle is wrapped to the interval between -180 + 180 degrees
%       ref.pha=wrapTo180(ref.pha);

% Second variant for vesions lower than 7.6

ref.pha=(180/pi)*angle(cos(ref.pha*pi/180)+1i*sin(ref.pha*pi/180));

ref.freq=waveform.inst_freq(ref.index);

% 1st Hypothesis
% ref.ROCOF=waveform.inst_ROCOF(ref.index);  

% 2nd Hypothesis
%  Npnts=max(size(ref.freq));
%  ref.ROCOF=([ref.freq(2)-ref.freq(1) diff(ref.freq)]+....
%  [diff(ref.freq) ref.freq(Npnts)-ref.freq(Npnts-1)])*param.Fs/2 ;

% 3rd  Hypothesis

% Application of a not casual digital filter

% OPTIONS:
% 1) Windows with a fixed number of elements (10)
NW=10; 
% For a variable NW (side elements of the window) tracked to 
% the frequency the  number can be computed by reference to the signal

% 2) Width of the window = reporting period
% NW=round(param.fs/param.Fs)-1; % Evaluation of the windows

% 3) Width of the window = period of the signal
% NW=round(param.fs/x.f0)-1; % Evaluation of the windows

%ref.ROCOF=filterF(NW,waveform.inst_ROCOF,ref.index);


