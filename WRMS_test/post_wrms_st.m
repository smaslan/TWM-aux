clc;
%clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);




% --- Post-processing of the WRMS algorithm uncertainty for single-tone inputs ---
% it uses N-dim data, reduces the dims to minimum needed by defining empirical relations...
% result is relatively small LUT and few equations

load('wrms_single_tone_bits.matsc','res','vr','p','s')

rr = [res{:}];

% get N-dim sizes:
par_n = vr.par_n;
par_n = par_n(par_n > 1);

% get error N-dim matrices:
dP = reshape([rr.dP],par_n);
dS = reshape([rr.dS],par_n);
dI = reshape([rr.dI],par_n);

% get no noise, max. bits errors:
%  note: limiting value - cannot be lower than that
dPn0 = repmat(dP(:,:,:,end,1),[1 1 1 par_n(end-1) par_n(end)]);
dSn0 = repmat(dS(:,:,:,end,1),[1 1 1 par_n(end-1) par_n(end)]);
dIn0 = repmat(dI(:,:,:,end,1),[1 1 1 par_n(end-1) par_n(end)]);

% get ref. value errors:
dPn1r1 = repmat(dP(1,:,:,end,end),[par_n(1) 1 1 par_n(end-1) par_n(end)]);
dSn1r1 = repmat(dS(1,:,:,end,end),[par_n(1) 1 1 par_n(end-1) par_n(end)]);
dIn1r1 = repmat(dI(1,:,:,end,end),[par_n(1) 1 1 par_n(end-1) par_n(end)]);

% load rms noise vector and cast it to its dimension:
rms_noise = reshape(p.rms_noise,[1 1 1 1 par_n(end)]);

% load bits vector and cast it to its dimension:
bits = reshape(p.b_bits,[1 1 1 par_n(end-1) 1]);

% load ab_rat vector and cast it to its dimension:
ab_rat = reshape(p.ab_rat,[par_n(1) 1 1 1 1]);

% effective rms noise comming from bit resolution:
%  note: empirically expanded
bits_noise_a = 1.5.*2.^-bits(end);
bits_noise_b = 1.5*ab_rat.*2.^-bits;

% total effective noise:
rms_noise_a = (bits_noise_a.^2 + rms_noise.^2).^0.5;
rms_noise_b = (bits_noise_b.^2 + rms_noise.^2).^0.5;

% combined relative dependence on rms noise and A amplitude:
rat_noise_a = (rms_noise_a./rms_noise(end))./(1./ab_rat(1));
% combined relative dependence on rms noise and B amplitude:
rat_noise_b = (rms_noise_b./rms_noise(end))./(ab_rat./ab_rat(1));

% estimated error components:
dP_est_a = dPn1r1.*rat_noise_a;
dP_est_b = dPn1r1.*rat_noise_b;
dS_est_a = dSn1r1.*rat_noise_a;
dS_est_b = dSn1r1.*rat_noise_b;
dI_est_b = dIn1r1.*rat_noise_b;

% total estimated errors:
dP_est = (dPn0.^2 + dP_est_a.^2 + dP_est_b.^2).^0.5;
dS_est = (dSn0.^2 + dS_est_a.^2 + dS_est_b.^2).^0.5;
dI_est = (dIn0.^2 + dI_est_b.^2).^0.5;

if false

    figure;
    semilogy(dP(:))
    hold on;
    semilogy(dP_est(:),'r')
    hold off; 
    
    figure;
    semilogy(dS(:))
    hold on;
    semilogy(dS_est(:),'r')
    hold off;
    
    figure;
    semilogy(dI(:))
    hold on;
    semilogy(dI_est(:),'r')
    hold off;

end


% extract data to be saved to the LUT:
dP_lut = dP(1,:,:,end,[1,par_n(end)])(:);
dS_lut = dS(1,:,:,end,[1,par_n(end)])(:);
dI_lut = dI(1,:,:,end,[1,par_n(end)])(:);

% build data vector:
eres = {struct()};
for c = 1:numel(dP_lut)    
    eres{c} = struct('dP',dP_lut(c),'dS',dS_lut(c),'dI',dI_lut(c));
end

% build list of parameter values (axes of dependence): 
ep.rms_noise = p.rms_noise([1,end]); % first and last noise value only
ep.f0_per = p.f0_per;
ep.fs_rat = p.fs_rat;

% -- build data vector descriptor
% axes of data vector:
evr.names = {'f0_per','fs_rat','rms_noise'};
% axes sizes:
evr.par_n = [numel(ep.f0_per) numel(ep.fs_rat) numel(ep.rms_noise)]; 


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
ax.rms_noise.min_lim = 'error';
ax.rms_noise.min_ovr = 0.99;
ax.rms_noise.max_lim = 'error';
ax.rms_noise.max_ovr = 1.05;
ax.rms_noise.scale = 'lin';
% define quantities interpolation behaviour:
qu = struct();
qu.dP.mult = 1.0; 
qu.dP.scale = 'log';
qu.dS.mult = 1.0; 
qu.dS.scale = 'log';
qu.dI.mult = 1.0; 
qu.dI.scale = 'log';

disp('Building LUT...');

lut = make_lut(eres,ep,evr,ax,qu);

% store reference a/b ratio:
lut.ref_ab_rat = ab_rat(1);
lut.info = 'Uncertainty LUT for windowed rms level and power algorithm. Made for blackmanharris() window.';

disp('Saving LUT...');

save('wrms_single_tone_unc.lut','-v7','lut')

%disp('Testing LUT...');
%interp_lut('test',lut,eres);




% --- Testing ---
% here the created LUT is used to calculate uncertainty and difference from source data is tested

dP_list = [];
for k = 1:1000

    % amplitdues of the channels:
    ab_rat_id = round(rand*(numel(ab_rat)-1)+1);
    amp_a = 1;
    amp_b = ab_rat(ab_rat_id);
    
    % rms noises of the channels:
    rms_noise_id = round(rand*(numel(rms_noise)-1)+1);
    rms_noise_a = rms_noise(rms_noise_id);
    rms_noise_b = rms_noise_a; % this will be different from the fist channel...
    
    % bits per pk-pk harmonic level:
    bits_id = round(rand*(numel(bits)-1)+1);
    bits_a = 32;
    bits_b = bits(bits_id);
    
    % sampling parameters:
    f0_per_id = round(rand*(numel(lut.ax.f0_per.values)-1)+1);
    f0_per = lut.ax.f0_per.values(f0_per_id);
    fs_rat_id = round(rand*(numel(lut.ax.fs_rat.values)-1)+1);
    fs_rat = lut.ax.fs_rat.values(fs_rat_id);
    
    % eval uncertainty:
    [dPx,dSx,sIx] = wrms_unc_st(lut, amp_a,amp_b, rms_noise_a,rms_noise_b, bits_a,bits_b, f0_per,fs_rat);
    
    % get reference value from source data:
    dPr = dP(ab_rat_id,f0_per_id,fs_rat_id,bits_id,rms_noise_id);
    dIr = dI(ab_rat_id,f0_per_id,fs_rat_id,bits_id,rms_noise_id);
    
    % store to list:
    dPx_list(k) = dPx;
    dPr_list(k) = dPr;

end

figure;
semilogy(dPr_list)
hold on;
semilogy(dPx_list,'r')
hold off;


