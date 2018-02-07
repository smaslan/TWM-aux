function [res] = thdwfft_sim_cell(sim)
% Part of non-coherent, windowed FFT, THD meter.
% 
% Wrapper for known THD simulator.
% This is called for each simulated set of parameters.
%

    % add alg. folder paths:
    cellfun(@addpath, sim.paths.paths, 'UniformOutput', false);
    
    
    % fetch correction data
    din = sim.corr;
    
    % merge corrections and setup
    din = struct([[fieldnames(din);fieldnames(sim.cfg)],[struct2cell(din);struct2cell(sim.cfg)]]'{:});
    
    % simulate THD waveforms
    [sig,fs,k1,a_ref] = thd_sim_wave(sim);
  
    % assign simulated data to the QWTB input
    din.Ts.v = 1/fs;
    din.y.v = sig;
        
    % backup paths (QWTB will remove algorithms path after calc.)
    path_backup = path();
    
    % QWTB calculation setup
    cset.verbose = sim.cfg.verbose;
    
    % execute algorithm once
    dout = qwtb('TWM-THDWFFT',din,cset);
    
    % restore paths
    path(path_backup);
    
    % return THD
    res.thd_ref = k1;
    res.thd = dout.thd.v;
    res.u_thd = dout.thd.u;
    res.thd_raw = dout.thd_raw.v;
    res.u_thd_raw = dout.thd_raw.u;
    
    % return amplitudes
    res.h_ref = a_ref;
    res.fh = dout.f.v;
    res.h = dout.h.v;
    res.u_h = dout.h.u;
    
    % return H2 amplitude
    res.h2_ref = a_ref(2);
    res.h2 = dout.h.v(2);
    res.u_h2 = dout.h.u(2);
    
    % return uncorrected H2 amplitude
    res.h2_raw = dout.h_raw.v(2);
    res.u_h2_raw = dout.h_raw.u(2);

end