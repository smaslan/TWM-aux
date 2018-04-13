clc;
clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);


    
    % --- multicore setup ---
    % multicore cores count
    mc_setup.cores = 100;
    % multicore method {'cellfun','parcellfun','multicore'}
    mc_setup.method = 'multicore';
    % multicore options: jobs grouping for 'parcellfun' 
    mc_setup.ChunksPerProc = 0;
    % multicore options: jobs grouping for 'multicore'
    mc_setup.min_chunk_size = 1;
    mc_setup.max_chunk_count = 50000;
    mc_setup.share_fld = [mfld filesep() 'mc_rubbish'];
    mc_setup.run_master_only = 0;
    mc_setup.run_slaves_only = 1;
    mc_setup.master_is_worker = 0;
    mc_setup.run_after_slaves = @coklbind2;
    
    % --- start slaves ---
    runmulticore(mc_setup.method,'',{},mc_setup.cores,mc_setup.share_fld,2,mc_setup);
        
    


