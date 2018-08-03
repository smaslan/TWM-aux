clc;
clear all;
close all;

% --- generate some parameters to variate:
% fundamental sine amplitude:
p.A0 = 1;
% fundamental sine periods in the window:
p.N0_per = linspace(3,5,50);
% sampling rate:
p.fs = 10000;
% fundamental frequency:
p.f0 = 50.3;
% second tone frequency:
p.fx = 150.1;
% second tone amplitude:
p.Ax = logspace(log10(100e-6),log10(1.0),7);
% monte-carlo cycles:
p.R = 50;
 

% initialize parameter variation lib for setup 'p'
[vr, par] = var_init(p);

% generate all combinations for of the input parameters 'p':
p_list = var_get_all_fast(par, vr, 5000, 1);
% now the 'p_list' contains vector of cells with all possible combinations
% of the vector parameters 'p'...


fprintf('Processing...\n')

% process the combinations in 'p_list' by function 'var_example_proc':
%  note: here should be prefferably some parallel processing method such as:
%   parfor - Matlab only
%   parcellfun - GNU Octave (Linux)
%   multicore - Matlab or GNU Octave
%   others?
% For simplicity we use single core processing 'cellfun':
res = cellfun(@var_example_proc,p_list,'UniformOutput',false);
% now 'res' should contain one result of the processing function
% for each combination in 'p_list'...


% --- 1) simple usage:

% now we can convert the results vector to dimensions mathing the variable parameter axes:
[r2d,ax_values,ax_names] = var_resize_result(res,vr,par);

% convert cell array of structs to single array 
max_err = reshape([[r2d{:}].e_rms],size(r2d));

% it is 2D dependence, so we can plot it as a mesh:
figure;
mesh(ax_values{2:-1:1},1e6*max(max_err,1e-8));
set(gca, 'xscale', 'log');
%set(gca, 'yscale', 'log');
%set(gca, 'zscale', 'log');
xlabel(str_insert_escapes(ax_names{2}));
ylabel(str_insert_escapes(ax_names{1}));
zlabel(str_insert_escapes('Windowed RMS error [ppm]'));



% --- 2) plotting dependencies along selected axes using var lib plots:

% plot dependence of 'res.e_rms' on the parameter 'par.N0_per' for parameter par.Ax(2)
var_plot_results_vect(vr,res,par,{'norm','e_rms','logy'},'N0_per','Ax',2);

% plot dependence of 'res.e_rms' on the parameter 'par.N0_per' for parameter par.Ax(3)
%  add limits plots defined as 'res.e_rms +- 2*res.s_rms'
var_plot_results_vect(vr,res,par,{'norm','e_rms','lim','s_rms',2,'logy'},'N0_per','Ax',3);

% plot dependence of 'res.e_rms' on the parameter 'par.Ax' for parameter par.N0_per(5)
var_plot_results_vect(vr,res,par,{'norm','e_rms','logx'},'Ax','N0_per',5);



