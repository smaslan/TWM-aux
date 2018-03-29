clc;
clear all;


M_list = [3:50];
j_list = logspace(-10,-5,20);
br_list = 6:24;
fsr_list = logspace(log10(10),log10(1000),30);

R = 500;

m_dfx = [];
s_dfx = [];
m_dAx = [];
s_dAx = [];
m_dpx = [];
s_dpx = [];

K = numel(br_list);
for k = 1:K

    fprintf('progress = %d%%\n',100*k/K);
    
    % maximum retries when non convergence or too high f deviation:
    max_try = 100;
                            
    % max. f deviation from estimate [-]:
    max_f_dev = 0.0005;
    
    % ADC bit resolution:
    %bits = 24;
    bits = br_list(k);
    
    % sampling time rms jitter [s]:
    jitt = 1e-9;
    %jitt = j_list(k);
    
    % harmonic frequencies:
    fh = 1e3;
    
    % periods of fundamental harmonic:
    M = 10;
    %M = M_list(k);
    
    % sampling rate ratio to fundamental:
    fs_rat = 100;
    %fs_rat = fsr_list(k);
    
    % sampling rate [Hz]:
    fs = fh*fs_rat;
    
    % samples count:
    N = round(fs/fh(1)*M);
    

    dhx = [];
    dox = [];
    dAx = [];
    dfx = [];
    for r = 1:R
        
        % generate time vector:
        t = [];
        t(:,1) = [0:N-1]/fs + jitt*randn(1,N);
        
        % generate random phases of harmonics:
        phi = rand(size(fh))*2*pi;
        ph0 = phi(1);
        
        % synthesize waveform:
        u = sin(2*pi*t.*fh + ph0);
        
        % add some random offset: 
        u = u + rand(1)*2^(bits-1);
        
        % round to ADC resolution:
        u = round(u*2^(bits-1))/2^(bits-1);
        
        for tr = 1:max_try
            % randomize initial guess:
            rand_t = randn(1)/fh(1)*(tr > 1);
            rand_f = max_f_dev*randn(1)*fh(1)*(tr > 1);
            
            % fit waveform:            
            [Ax, fx, phx, ox] = FPNLSF(t+rand_t,u,fh(1)*(1+rand_f),0);
            phx = phx - fx*rand_t*2*pi;
                        
            if ~isinf(fx) && abs(fx/fh(1)-1) < max_f_dev
                % result possibly ok - leave 
                break;
            elseif tr == max_try
                disp('No convergence even after all retires! Dunno what to do now...');
                plot(u)
%             else
%                 abs(fx/fh(1)-1)
            end
            % retry because we got no convergence or too high phase deviation
                     
        end
        
        
        % calculate phase deviation (wrapped to +-pi):
        dhx = (phx - ph0);
        dhx(r) = mod(dhx+pi,2*pi)-pi;
        
        % deviation of ampl.:
        dAx(r) = Ax - 1;
        
        % deviation of freq.:
        dfx(r) = fx/fh(1)-1;
        
        % deviation of offset.:
        dox(r) = ox - 0;
    
    end
    
    fid = 1:numel(dhx);
    
    m_dfx(k) = mean(dfx(fid));
    s_dfx(k) = std(dfx(fid));
    
    m_dAx(k) = mean(dAx(fid));
    s_dAx(k) = std(dAx(fid));
    
    m_dpx(k) = mean(dhx(fid));
    s_dpx(k) = std(dhx(fid));

end

 






