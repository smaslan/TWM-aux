clc;
clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

% -----------------------------------------------------------------------------------------------------------------------
% Sensitivity analysis for JV's FFT filter 'sFreqDep_PG_Comp.m'
% -----------------------------------------------------------------------------------------------------------------------

% samples count (2^N):
%  note: 10 means it will simulate random sample counts from (2^9)+1 to 2^10 
s.N_pow = [10:2:20];
% harmonic frequency range relative to fs:
s.f0_rat_min_bin = 6;
s.f0_rat_max = [0.02 0.05 0.1 0.2 0.3 0.4 0.45];
% harmonic amplitude:
s.f0_amp = 1;
% randomize phase range [+- rad]:
s.f0_phi_rnd = pi;
% filter max phase [+- rad]:
s.ff_max_phi = [0 logspace(-6,-0,10)]*pi;
% filter max real amp dev [+- V/V]:
s.ff_max_amp = [0 logspace(-6,log10(0.5),10)]*pi;
% bits count per f0_amp:
s.bits = [24 32];
% rms noise:
s.rms_noise = [0 logspace(-7,-3,5)];
% tests per setting:
s.cycles = 10000;


% --- initialize simulator/generate jobs ---
[vr,p] = var_init(s);
p_list = var_get_all_fast(p,vr,5000,1);



% --- multicore setup ---
% multicore cores count
mc_setup.cores = 576/2;
% multicore method {'cellfun','parcellfun','multicore'}
mc_setup.method = 'multicore';
% multicore options: jobs grouping for 'parcellfun' 
mc_setup.ChunksPerProc = 0;
% multicore options: jobs grouping for 'multicore'
mc_setup.share_fld = [mfld filesep() 'mc_rubbish'];
mc_setup.min_chunk_size = 1;
if ispc
    % windoze - small CPU
    mc_setup.master_is_worker = 1;    
    mc_setup.max_chunk_count = 200;
    mc_setup.run_master_only = 1;    
else
    % cokl supercomputer - large CPU
    mc_setup.master_is_worker = 0;    
    mc_setup.max_chunk_count = 20000;
    mc_setup.run_master_only = 0;
    mc_setup.run_after_slaves = @coklbind2;
end

% --- process batch ---
res = runmulticore(mc_setup.method,@proc_FFTF,p_list,mc_setup.cores,mc_setup.share_fld,2,mc_setup);
% update result count in jobs list
vr.res_n = length(res);


save('fftf_test.matsc','-v7','res','vr','p','s')



