clc;
%clear all;

mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);





% --- WRMS algorithm uncertainty for spurr effects ---
% it uses N-dim data, reduces the dims to minimum needed by defining empirical relations...
% result is relatively small LUT and few equations

%load('wrms_spurr_test.matsc','res','vr','p','s')

for k = 1:9

    f0_per_id = k;
    fs_rat_id = 5;
    
    [r,v,s_freq] = var_get_results_vect(res,p,vr,'s_freq','f0_per',f0_per_id,'fs_rat',fs_rat_id);
    
    f0 = 1;
    f0_per = p.f0_per(f0_per_id);
    fs = p.fs_rat(fs_rat_id);
    
    N = f0_per*fs
    
    
    spurr_f = s_freq.*(fs/2 - f0) + f0;
    
    dft_step = fs/N;
    
    dft_space = spurr_f/dft_step;
    
    
    loglog(dft_space,v.dP)
    hold on;

end
hold off;




