clc;
%clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);




% --- Post-processing of the WRMS algorithm uncertainty for spurr effects ---
% it uses N-dim data, reduces the dims to minimum needed by defining empirical relations...
% result is relatively small LUT and few equations

%load('wrms_spurr_test.matsc','res','vr','p','s')

rr = [res{:}];

% get N-dim sizes:
par_n = vr.par_n;
par_n = par_n(par_n > 1);

% get error N-dim matrices:
dP = reshape([rr.dP],par_n);
dS = reshape([rr.dS],par_n);
dI = reshape([rr.dI],par_n);

% get minimum error case values:
%  note: limiting value - cannot be lower than that
dPn0 = repmat(dP(end,:,:,1,:),[par_n(1) 1 1 par_n(end-1) 1]);
dSn0 = repmat(dS(end,:,:,1,:),[par_n(1) 1 1 par_n(end-1) 1]);
dIn0 = repmat(dI(end,:,:,1,:),[par_n(1) 1 1 par_n(end-1) 1]);

% get ref. value errors:
%dIn1r1 = repmat(dI(1,:,:,end,:),[par_n(1) 1 1 par_n(end-1) 1]);
%dIn1r1 = repmat(dI(:,:,:,end,:),[1 1 1 par_n(end-1) 1]);
dIn1r1 = repmat(dI(end,:,:,:,:),[par_n(1) 1 1 1 1]);
dPn1r1 = repmat(dP(end,:,:,:,:),[par_n(1) 1 1 1 1]);
dSn1r1 = repmat(dS(end,:,:,:,:),[par_n(1) 1 1 1 1]);

% load ab_rat vector and cast it to its dimension:
ab_rat = reshape(p.ab_rat,[par_n(1) 1 1 1 1]);

% load spurr amplitude vector:
s_amp = reshape(p.s_amp,[1 1 1 par_n(end-1) 1]);

% combined relative dependence on rms noise and B amplitude:
%rat_gain_b = (s_amp./s_amp(end))./(ab_rat./ab_rat(1));
%rat_gain_b = (s_amp./s_amp(end));
rat_gain_b = 1./(ab_rat./ab_rat(end));
rat_gain_a = 1./(1./ab_rat(end));

% ratio contributions:
dP_est_a = dPn1r1.*rat_gain_a;
dP_est_b = dPn1r1.*rat_gain_b;
dS_est_a = dSn1r1.*rat_gain_a;
dS_est_b = dSn1r1.*rat_gain_b;
dI_est_b = dIn1r1.*rat_gain_b;

% total estimated errors:
dP_est = (dPn0.^2 + dP_est_b.^2 + dP_est_a.^2).^0.5;
dS_est = (dPn0.^2 + dS_est_b.^2 + dS_est_a.^2).^0.5;
dI_est = (dIn0.^2 + dI_est_b.^2).^0.5;


if true

    figure;
    semilogy(dP(:))
    hold on;
    semilogy(dP_est(:),'r')
    hold off;
    title('dP'); 
    
    figure;
    semilogy(dP_est(:)./dP(:))
    title('dP ratio');
    
    figure;
    semilogy(dS(:))
    hold on;
    semilogy(dS_est(:),'r')
    hold off;
    title('dS');
    
    figure;
    semilogy(dI(:))
    hold on;
    semilogy(dI_est(:),'r')
    hold off;
    title('dI')
    
    figure;
    semilogy(dI_est(:)./dI(:))
    title('dI ratio');

end

return


% --- generate LUT A ---

% extract data to be saved to the LUT:
dP_lut = dP(1,:,:,:,:)(:);
dS_lut = dS(1,:,:,:,:)(:);
dI_lut = dI(1,:,:,:,:)(:);

% build data vector:
eres = {struct()};
for c = 1:numel(dP_lut)    
    eres{c} = struct('dP',dP_lut(c),'dS',dS_lut(c),'dI',dI_lut(c));
end

% build list of parameter values (axes of dependence): 
ep = struct();
ep.f0_per = p.f0_per;
ep.fs_rat = p.fs_rat;
ep.s_amp = p.s_amp;
ep.s_freq = p.s_freq;

% -- build data vector descriptor
% axes of data vector:
evr.names = {'f0_per','fs_rat','s_amp','s_freq'};
% axes sizes:
evr.par_n = [numel(ep.f0_per) numel(ep.fs_rat) numel(ep.s_amp) numel(ep.s_freq)]; 


% define axes limits and behaviour:
ax = struct();
ax.f0_per.min_lim = 'error';
ax.f0_per.min_ovr = 0.99;
ax.f0_per.max_lim = 'const';
ax.f0_per.max_ovr = 1.05;
ax.f0_per.scale = 'log';
ax.fs_rat.min_lim = 'const';
ax.fs_rat.min_ovr = 0.99;
ax.fs_rat.max_lim = 'const';
ax.fs_rat.max_ovr = 1.05;
ax.fs_rat.scale = 'log';
ax.s_amp.min_lim = 'const';
ax.s_amp.min_ovr = 0.99;
ax.s_amp.max_lim = 'const';
ax.s_amp.max_ovr = 1.05;
ax.s_amp.scale = 'log';
ax.s_freq.min_lim = 'const';
ax.s_freq.min_ovr = 0.99;
ax.s_freq.max_lim = 'const';
ax.s_freq.max_ovr = 1.05;
ax.s_freq.scale = 'log';

% define quantities interpolation behaviour:
qu = struct();
qu.dP.mult = 1.0; 
qu.dP.scale = 'log';
qu.dS.mult = 1.0; 
qu.dS.scale = 'log';
qu.dI.mult = 1.0; 
qu.dI.scale = 'log';

disp('Building LUT A...');

lut.lut_ref = make_lut(eres,ep,evr,ax,qu);

% store reference a/b ratio:
lut.ref_ab_rat = ab_rat(1);
lut.info = 'Uncertainty LUT for windowed rms level and power algorithm - spurr effect. Made for blackmanharris() window.';


% --- generate LUT B ---

% extract data to be saved to the LUT:
dP_lut = dP(end,:,:,1,:)(:);
dS_lut = dS(end,:,:,1,:)(:);
dI_lut = dI(end,:,:,1,:)(:);

% build data vector:
eres = {struct()};
for c = 1:numel(dP_lut)    
    eres{c} = struct('dP',dP_lut(c),'dS',dS_lut(c),'dI',dI_lut(c));
end

% build list of parameter values (axes of dependence): 
ep = struct();
ep.f0_per = p.f0_per;
ep.fs_rat = p.fs_rat;
ep.s_freq = p.s_freq;

% -- build data vector descriptor
% axes of data vector:
evr.names = {'f0_per','fs_rat','s_freq'};
% axes sizes:
evr.par_n = [numel(ep.f0_per) numel(ep.fs_rat) numel(ep.s_freq)]; 


% define axes limits and behaviour:
ax = struct();
ax.f0_per.min_lim = 'error';
ax.f0_per.min_ovr = 0.99;
ax.f0_per.max_lim = 'const';
ax.f0_per.max_ovr = 1.05;
ax.f0_per.scale = 'log';
ax.fs_rat.min_lim = 'const';
ax.fs_rat.min_ovr = 0.99;
ax.fs_rat.max_lim = 'const';
ax.fs_rat.max_ovr = 1.05;
ax.fs_rat.scale = 'log';
ax.s_freq.min_lim = 'const';
ax.s_freq.min_ovr = 0.99;
ax.s_freq.max_lim = 'const';
ax.s_freq.max_ovr = 1.05;
ax.s_freq.scale = 'log';

% define quantities interpolation behaviour:
qu = struct();
qu.dP.mult = 1.0; 
qu.dP.scale = 'log';
qu.dS.mult = 1.0; 
qu.dS.scale = 'log';
qu.dI.mult = 1.0; 
qu.dI.scale = 'log';

disp('Building LUT B...');

lut.lut_min = make_lut(eres,ep,evr,ax,qu);


disp('Saving LUT...');

save('wrms_spurr_unc.lut','-v7','lut')

%disp('Testing LUT...');
%interp_lut('test',lut,eres);




% --- Testing ---
% here the created LUT is used to calculate uncertainty and difference from source data is tested

disp('Testing LUT...');

tic
dPr_list = [];
dPx_list = [];
dIr_list = [];
dIx_list = [];
for k = 1:100

    % amplitdues of the channels:
    ab_rat_id = round(rand*(numel(ab_rat) - 1) + 1);
    amp_a = 1;
    amp_b = ab_rat(ab_rat_id);
    
    % rms noises of the channels:
    f_spurr_id = round(rand*(numel(p.s_freq) - 1) + 1);
    f_spurr = p.s_freq(f_spurr_id);
        
    % bits per pk-pk harmonic level:
    a_spurr_id = round(rand*(numel(p.s_amp) - 1) + 1);
    a_spurr = p.s_amp(a_spurr_id);
    
    % sampling parameters:
    f0_per_id = round(rand*(numel(lut.ax.f0_per.values)-1)+1);
    f0_per = lut.ax.f0_per.values(f0_per_id);
    fs_rat_id = round(rand*(numel(lut.ax.fs_rat.values)-1)+1);
    fs_rat = lut.ax.fs_rat.values(fs_rat_id);
    
    % eval uncertainty:
    [dPx,dSx,dIx] = wrms_unc_spurr(lut, amp_a,amp_b, f_spurr,a_spurr, f0_per,fs_rat);
    
    % get reference value from source data:
    dPr = dP(ab_rat_id,f0_per_id,fs_rat_id,a_spurr_id,f_spurr_id);
    dIr = dI(ab_rat_id,f0_per_id,fs_rat_id,a_spurr_id,f_spurr_id);
    
    % store to list:
    dPx_list(k) = dPx;
    dPr_list(k) = dPr;
    dIx_list(k) = dIx;
    dIr_list(k) = dIr;

end
toc

figure;
semilogy(dPr_list)
hold on;
semilogy(dPx_list,'r')
hold off;
title('P');

figure;
semilogy(dIr_list)
hold on;
semilogy(dIx_list,'r')
hold off;
title('I');



