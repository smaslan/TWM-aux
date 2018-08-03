function [res] = var_example_proc(par)
% This is example function for 'var' lib demonstration.
% It calculates maximum error of windowed RMS algorithm for given setup:
%  par.fs     = sampling rate
%  par.A0     = amplitude of fundamental component  
%  par.f0     = frequency of fundamental component
%  par.N0_per = periods of fundamental component (window size)
%  par.fx     = frequency of second tone
%  par.Ax     = amplitude of second tone
%  par.R      = monte-carlo tries
%
% This function will try 'R' times to repeat the calculation with given setup
% and random phase of the signal. Then it will calculate maximum error of all
% the tries.
%
% License:
% --------
% This is part of VAR library for automatic multidim. variation of simulation parameters.
% (c) 2018, Stanislav Maslan, s.maslan@seznam.cz
% The script is distributed under MIT license, https://opensource.org/licenses/MIT 

    % samples count per waveform:
    N = round(par.N0_per/par.f0*par.fs);
    
    % generate time vector:
    t(:,1) = [0:N-1]/par.fs*2*pi;
    
    % generate Hann window:
    w(:,1) = 0.5 - 0.5*cos([0:N-1]/N*2*pi);
        
    % window RMS value:
    w_rms = mean(w.^2)^0.5;
    
    
    % generate sine wave with random phase (generate all monte-carlo cycles at once):
    %  note: this is crippled to be able to run in Matlab < 2016b, was: 
    %        u = par.A*sin(t*par.f + rand(1,par.R)*2*pi);
    u = par.A0*sin(bsxfun(@plus,t*par.f0,rand(1,par.R)*2*pi));
    % generate second tone
    u = u + par.Ax*sin(bsxfun(@plus,t*par.fx,rand(1,par.R)*2*pi));
    
    % apply window to it:         
    %  u = u.*w;
    u = bsxfun(@times,u,w);
    
    % calculate RMS (one per monte carlo cycle):
    u_rms = mean(u.^2,1).^0.5/w_rms;
    
    % reference RMS value:
    r_rms = sum(0.5*[par.A0,par.Ax].^2)^0.5;
    
    % store maximum RMS error:
    res.e_rms = max(abs(u_rms - r_rms)./r_rms);   
    
    % store stdev of RMS error:
    res.s_rms = std(abs(u_rms - r_rms)./r_rms);

end

