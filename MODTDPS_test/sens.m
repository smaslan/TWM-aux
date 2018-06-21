clc;
clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

warning('off');

id = 0;

id++;
cfg{id}.name = 'Sine - corrections';
cfg{id}.file = 'sine_corr_v2.matsc';
cfg{id}.is_rect = 0;
cfg{id}.comp_err = 1;

% id++;
% cfg{id}.name = 'Sine - no corrections';
% cfg{id}.file = 'sine_ncorr.matsc';
% cfg{id}.is_rect = 0;
% cfg{id}.comp_err = 0;
% 
% id++;
% cfg{id}.name = 'Rectangular - no corrections';
% cfg{id}.file = 'rect_ncorr.matsc';
% cfg{id}.is_rect = 1;
% cfg{id}.comp_err = 0;


for k = 1:numel(cfg)

    clear res;
    clear p;


    disp(['Running: ' cfg{k}.name]);

    s = struct();
    % ADC bit resolution of the pk-pk amplitude:
    s.bits = unique(round(logspace(log10(4),log10(24),6)));
    % sampling time rms jitter [s]:
    s.jitt = logspace(-9,-2,5);
    % modulate/carrier freq ratio:
    s.fmf0_rat = logspace(log10(0.01),log10(0.33),8);
    % modulating depth:
    s.modd = logspace(log10(0.01),log10(0.99),8);
    % sampling rate to carrier freq ratio:
    s.fsf0_rat = logspace(log10(10),log10(100),5);
    % modulating wave shape:
    s.is_rect = cfg{k}.is_rect;
    % self-compensate error:
    s.comp_err = cfg{k}.comp_err;
    % periods count of the modulating signal:
    s.fm_per = logspace(log10(3),log10(30),6);
    % offset randomizing:
    s.ofs = 0.01;
    % sfdr values:
    s.sfdr = [120 80 60 30];
    % sfdr non harmonic:
    s.sfdr_nc = 1;
     
    % mcc:
    s.R = 1000;
    s.R_max = 10000;
        
    
    % --- initialize simulator/generate jobs ---
    [vr,p] = var_init(s);
    p_list = var_get_all_fast(p,vr,5000,1);
    
    fprintf('Count: %d\n',numel(p_list));
    
    % --- multicore setup ---
    % multicore cores count
    mc_setup.cores = 300;
    % multicore method {'cellfun','parcellfun','multicore'}
    mc_setup.method = 'multicore';
    % multicore options: jobs grouping for 'parcellfun' 
    mc_setup.ChunksPerProc = 0;
    % multicore options: jobs grouping for 'multicore'
    mc_setup.min_chunk_size = 1;
    mc_setup.max_chunk_count = 50000;
    mc_setup.share_fld = [mfld filesep() 'mc_rubbish'];
    mc_setup.run_master_only = 0;
    mc_setup.master_is_worker = 0;
    mc_setup.run_after_slaves = @coklbind2;
    
    % --- process batch ---
    res = runmulticore(mc_setup.method,@proc_MODTDPS,p_list,mc_setup.cores,mc_setup.share_fld,2,mc_setup);
    % update result count in jobs list
    vr.res_n = length(res);
    
    
    save(cfg{k}.file,'-v6','res','vr','p','s');

end


