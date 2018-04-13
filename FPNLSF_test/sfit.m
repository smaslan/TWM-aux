
mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

disp('Loading data...');

% load results:
load('test.matsc','-mat');

disp('Preprocessing data...');

% merge deviations and uncertainties:
rr = {struct()};
for k = 1:numel(res)   
   rr{k}.dpx = res{k}.s_dpx;
   rr{k}.dfx = res{k}.s_dfx;
   rr{k}.dAx = res{k}.s_dAx;
   rr{k}.dox = res{k}.s_dox;
end

disp('Building LUT...');

% define axes limits and behaviour:
ax = struct();
ax.bits.min_lim = 'error';
ax.bits.min_ovr = 0.99;
ax.bits.max_lim = 'const';
ax.bits.max_ovr = 1.05;
ax.bits.scale = 'lin';
ax.f0_per.min_lim = 'error';
ax.f0_per.min_ovr = 0.99;
ax.f0_per.max_lim = 'const';
ax.f0_per.max_ovr = 1.05;
ax.f0_per.scale = 'log';
ax.fs_rat.min_lim = 'error';
ax.fs_rat.min_ovr = 0.99;
ax.fs_rat.max_lim = 'const';
ax.fs_rat.max_ovr = 1.05;
ax.fs_rat.scale = 'log';
ax.jitt.min_lim = 'const';
ax.jitt.min_ovr = 0.99;
ax.jitt.max_lim = 'error';
ax.jitt.max_ovr = 1.05;
ax.jitt.scale = 'log';
ax.sfdr.min_lim = 'const';
ax.sfdr.min_ovr = 0.99;
ax.sfdr.max_lim = 'error';
ax.sfdr.max_ovr = 1.05;
ax.sfdr.scale = 'lin';
% define quantities interpolation behaviour:
qu = struct();
qu.dpx.mult = 1.0; 
qu.dpx.scale = 'log';
qu.dfx.mult = 1.0;
qu.dfx.scale = 'log';
qu.dAx.mult = 1.0;
qu.dAx.scale = 'log';
qu.dox.mult = 1.0;
qu.dox.scale = 'log';


lut = make_lut(rr,p,vr,ax,qu);

disp('Saving LUT...');

save('unc.lut','-v7','lut')

disp('Testing LUT...');

interp_lut('test',lut,rr);




