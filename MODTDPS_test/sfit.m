
mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

disp('Loading data...');

% load results:
%load('sine_corr_min.matsc','-mat');

disp('Preprocessing data...');

% merge deviations and uncertainties:
% rr = {struct()};
% for k = 1:numel(res)   
%    %rr{k}.dofs = res{k}.s_dofs;
%    rr{k}.df0 = res{k}.s_df0;
%    rr{k}.dfm = res{k}.s_dfm;
%    rr{k}.dA0 = res{k}.s_dA0;
%    rr{k}.dmod = res{k}.s_dmod;   
% end

disp('Building LUT...');

% define axes limits and behaviour:
ax = struct();
ax.bits.min_lim = 'error';
ax.bits.min_ovr = 0.99;
ax.bits.max_lim = 'const';
ax.bits.max_ovr = 1.05;
ax.bits.scale = 'lin';

ax.jitt.min_lim = 'const';
ax.jitt.min_ovr = 0.99;
ax.jitt.max_lim = 'error';
ax.jitt.max_ovr = 1.05;
ax.jitt.scale = 'log';

ax.fmf0_rat.min_lim = 'const';
ax.fmf0_rat.min_ovr = 0.99;
ax.fmf0_rat.max_lim = 'error';
ax.fmf0_rat.max_ovr = 1.02;
ax.fmf0_rat.scale = 'log';

ax.fsf0_rat.min_lim = 'error';
ax.fsf0_rat.min_ovr = 0.99;
ax.fsf0_rat.max_lim = 'const';
ax.fsf0_rat.max_ovr = 1.05;
ax.fsf0_rat.scale = 'log';

ax.modd.min_lim = 'error';
ax.modd.min_ovr = 0.999;
ax.modd.max_lim = 'error';
ax.modd.max_ovr = 1.001;
ax.modd.scale = 'log';

ax.fm_per.min_lim = 'error';
ax.fm_per.min_ovr = 0.99;
ax.fm_per.max_lim = 'const';
ax.fm_per.max_ovr = 1.05;
ax.fm_per.scale = 'log';

ax.sfdr.min_lim = 'const';
ax.sfdr.min_ovr = 0.99;
ax.sfdr.max_lim = 'error';
ax.sfdr.max_ovr = 1.05;
ax.sfdr.scale = 'lin';



% define quantities interpolation behaviour:
qu = struct();
%qu.dofs.mult = 1.5; 
%qu.dofs.scale = 'log';
%qu.dofs.vmax = inf;
qu.df0.mult = 1.5;
qu.df0.scale = 'log';
qu.df0.vmax = 1.0;
qu.dfm.mult = 1.5;
qu.dfm.scale = 'log';
qu.dfm.vmax = 1.0;
qu.dA0.mult = 1.5;
qu.dA0.scale = 'log';
qu.dA0.vmax = 1.0;
qu.dmod.mult = 1.5;
qu.dmod.scale = 'log';
qu.dmod.vmax = 1.0;

lut = make_lut(rr,p,vr,ax,qu);

disp('Saving LUT...');

save('sine_corr_unc.lut','-v7','lut')

disp('Testing LUT...');

interp_lut('test',lut,rr);




