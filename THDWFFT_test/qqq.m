clc;
clear all;
close all;

% this folder
mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% TWM/octprog/
ofld = [mfld filesep() '..' filesep() '..' filesep() 'TWM' filesep() 'octprog'];

% add TWM/octprog
addpath(ofld);

% add QWTB and INFO-STRINGS paths
%addpath([ofld filesep() 'info']);
qwtb_path = [ofld filesep() 'qwtb'];
addpath(qwtb_path);

% add path to the tested algorithm
thd_path = [ofld filesep() 'qwtb' filesep() 'alg_TWM-THDWFFT'];
addpath(thd_path);



% --- correction data ---
cal.adc_gain.v = [1.0;0.5];
cal.adc_gain.u = 0*[0.01;0.01];
cal.adc_gain_f.v = [0;1e5];
cal.adc_gain_a.v = [];

cal.tr_gain.v = [1.0];
cal.tr_gain.u = [0.0];
cal.tr_gain_f.v = [];
cal.tr_gain_a.v = [];

cal.adc_sfdr.v = [180];
cal.adc_sfdr_f.v = [];
cal.adc_sfdr_a.v = [];

cal.tr_sfdr.v = [180];
cal.tr_sfdr_f.v = [];
cal.tr_sfdr_a.v = [];

cal.adc_nrng.v = 1.0;
cal.adc_bits.v = 24;

% Restore orientations of the input vectors to originals (before passing via QWTB)
cal.y.v = ones(10,1); % fake data vector just to make following function work!
[cal,scfg] = qwtb_restore_twm_input_dims(cal,1);
% Rebuild TWM style correction tables (just for more convenient calculations):
tab = qwtb_restore_correction_tables(cal,scfg);



% --- algorithm setup ---
% plot spectrum?
cfg.plot.v = 0;
% harmonics count to analyze
cfg.H.v = 10;
% verbose
cfg.verbose.v = 0;
% initial guess of the fundamental frequency (comment if autodetect needed)
%cfg.f0.v = 1e3;
% fix scalloping error?
cfg.scallop_fix.v = 1;



% --- simulator setup ---
% repetitions
sim.reps = 0;[1:10];
% sampling rate [Hz]
sim.fs = 50e3;
% fundamental freq [Hz]
sim.f0 = 1.0e3;
%sim.f0 = cfg.f0.v;

% samples count
sim.sample_count = round(sim.fs/sim.f0*200);
%sim.sample_count = round(sim.fs/sim.f0*logspace(log10(50),log10(1000),50));

% fundamental freq [Hz] - sweep
f_bin_step = sim.fs/sim.sample_count;
sim.f0 = sim.f0 + 0.6*linspace(-f_bin_step,f_bin_step,51);

% repeated measurements (averages) count
sim.avg_count = 3;
% fundamental amplitude [V]
sim.A0 = 0.9;
% --- to generate:
    % 1) desired THD (fundamental referenced) 
    %sim.k1 = 0.005; %logspace(log10(0.0001),log10(10),50);
    % 2) or fixed, identical amplitudes [V]
    sim.A = 0.1; %logspace(log10(1e-7),log10(1e-5),20);
    % 3) or random amplitudes in logspace, range from-to [V]
    %sim.A_min = 1e-6;
    %sim.A_max = 1000e-6;
% harmonics count to generate (including fundamental)
sim.H = cfg.H.v;
% ADC noise level in DFT spectrum [V]
sim.adc_noise_lev = 1e-6;logspace(log10(1e-7),log10(10e-6),0);
% enable randomization of quantities with uncertainties
sim.randomize = 0;
% copy algorithm calculation setup
sim.corr = cal;
sim.tab = tab;
sim.cfg = cfg;
% paths to THD alg. components (used for multicore processing)
sim.paths.paths = {qwtb_path; thd_path};

% --- initialize simulator/generate jobs ---
[vr,p] = var_init(sim);
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
mc_setup.cores = 4;
% multicore method {'cellfun','parcellfun','multicore'}
mc_setup.method = 'multicore';
% multicore options: jobs grouping for 'parcellfun' 
mc_setup.ChunksPerProc = 0;
% multicore options: jobs grouping for 'multicore'
mc_setup.min_chunk_size = 1;
mc_setup.max_chunk_count = 100;
mc_setup.share_fld = [mfld filesep() 'mc_rubbish'];
mc_setup.run_master_only = 1;

% --- process batch ---
res = runmulticore(mc_setup.method,@thdwfft_sim_cell,p_list,mc_setup.cores,mc_setup.share_fld,2,mc_setup);
% update result count in jobs list
vr.res_n = length(res);



