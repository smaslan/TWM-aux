clc;
clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);


% I/U amplitude ratio:
s.ab_rat = logspace(log10(0.01),log10(1.0),5);
s.ab_rat_rnd = 1e-3;
% I-U pahse shift:
s.ab_phi = 0;
s.ab_phi_rnd = pi;
% periods of fundamental harmonic:
s.f0_per = [3,4,5,6,7,9,12,20,50,100];
% sampling rate ratio to f0 ratio:
s.fs_rat = logspace(log10(7),log10(100),10);
% spurr amplitude:
s.s_amp = logspace(log10(0.001),log10(1.0),9);
% spurr freq to fs/2 ratio:
s.s_freq = logspace(log10(0.001),log10(0.9),15);
% spurr analysis mode?
s.spurr_mode = 1;
% bits count:
s.a_bits = 32;
s.b_bits = 32;logspace(log10(4),log10(32),8);
% rms noise:
s.rms_noise = 0;[0 logspace(-7,-3,5)];
% tests per setting:
s.cycles = 300;





% --- initialize simulator/generate jobs ---
[vr,p] = var_init(s);
p_list = var_get_all_fast(p,vr,5000,1);

% --- multicore setup ---
% multicore cores count
mc_setup.cores = 300;
% multicore method {'cellfun','parcellfun','multicore'}
mc_setup.method = 'multicore';
% multicore options: jobs grouping for 'parcellfun' 
mc_setup.ChunksPerProc = 0;
% multicore options: jobs grouping for 'multicore'
mc_setup.share_fld = [mfld filesep() 'mc_rubbish'];
mc_setup.min_chunk_size = 1;
mc_setup.master_is_worker = 0;
if ispc
    % windoze - small CPU    
    mc_setup.max_chunk_count = 1000;
    mc_setup.run_master_only = 1;    
else
    % cokl supercomputer - large CPU    
    mc_setup.max_chunk_count = 10000;
    mc_setup.run_master_only = 0;
    mc_setup.run_after_slaves = @coklbind2;
end

% --- process batch ---
res = runmulticore(mc_setup.method,@proc_WRMS,p_list,mc_setup.cores,mc_setup.share_fld,2,mc_setup);
% update result count in jobs list
vr.res_n = length(res);


save('wrms_spurr_test.matsc','-v7','res','vr','p','s')



