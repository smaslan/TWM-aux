clc;
clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);


% I/U amplitude ratio:
s.ab_rat = logspace(log10(0.01),log10(1.0),15);
% I-U pahse shift:
s.ab_phi = 0;[linspace(0,0.5,20) linspace(0.51,0.99,20)]*pi;
s.ab_phi_rnd = pi;
% periods of fundamental harmonic:
s.f0_per = [3,4,5,6,7,9,11,15,20,50,100];
% sampling rate ratio to f0 ratio:
s.fs_rat = logspace(log10(7),log10(100),12);
% spurr amplitude:
s.s_amp = 0;
% spurr freq to fundamental freq:
s.s_freq = 0;
% bits count:
s.a_bits = 32;
s.b_bits = logspace(log10(4),log10(32),9);
% rms noise:
s.rms_noise = [0 logspace(-7,-3,5)];
% tests per setting:
s.cycles = 50000;





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
    mc_setup.max_chunk_count = 200;
    mc_setup.run_master_only = 1;    
else
    % cokl supercomputer - large CPU    
    mc_setup.max_chunk_count = 20000;
    mc_setup.run_master_only = 0;
    mc_setup.run_after_slaves = @coklbind2;
end

% --- process batch ---
res = runmulticore(mc_setup.method,@proc_WRMS,p_list,mc_setup.cores,mc_setup.share_fld,2,mc_setup);
% update result count in jobs list
vr.res_n = length(res);


save('wrms_single_tone_bits.matsc','-v7','res','vr','p','s')



