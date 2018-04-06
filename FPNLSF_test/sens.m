clc;
clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);



% maximum retries when non convergence or too high f deviation:
s.cfg.max_try = 100;
% max. f deviation from estimate [-]:
s.cfg.max_f_dev = 500e-6;
% ADC bit resolution:
s.bits = unique(round(logspace(log10(4),log10(24),8)));
% sampling time rms jitter [s]:
s.jitt = logspace(-9,-2,9);
% periods of fundamental harmonic:
s.f0_per = [10,20,50,100];
% sampling rate ratio to f0 ratio:
s.fs_rat = logspace(log10(10),log10(1000),10);
% sfdr values:
s.sfdr = [180 120 80 40 30];
% sfdr non harmonic:
s.sfdr_nc = 1;
% mcc:
s.R = 1000;



% --- initialize simulator/generate jobs ---
[vr,p] = var_init(s);
if vr.var_n == 1
    % enable spectrum logging if only one simulation enabled
    s.save_spec = 1;
    p.s = s;
    % reinit simulator structure
    [vr,p] = var_init(p);
endif
p_list = var_get_all(p,vr,5000,1);

% --- multicore setup ---
% multicore cores count
mc_setup.cores = 150;
% multicore method {'cellfun','parcellfun','multicore'}
mc_setup.method = 'multicore';
% multicore options: jobs grouping for 'parcellfun' 
mc_setup.ChunksPerProc = 0;
% multicore options: jobs grouping for 'multicore'
mc_setup.min_chunk_size = 1;
mc_setup.max_chunk_count = 2000;
mc_setup.share_fld = [mfld filesep() 'mc_rubbish'];
mc_setup.run_master_only = 0;
mc_setup.run_after_slaves = @coklbind2;

% --- process batch ---
res = runmulticore(mc_setup.method,@proc_FPNLSF,p_list,mc_setup.cores,mc_setup.share_fld,2,mc_setup);
% update result count in jobs list
vr.res_n = length(res);


save('cokl_long.mat','res','vr','p','s')



