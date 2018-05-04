function [r] = proc_WRMS(p)

    warning('off');

    
    dA = [];
    dP = [];
    dN = [];
    dR = [];
    for k = 1:p.cycles
    
        % random size from 2^(p.N_pow-1)+1 to 2^(p.N_pow) 
        N = 2^p.N_pow - round((2^(p.N_pow-1) - 1)*rand(1));
        N = floor(N/2)*2;
        
        % size of the FFT filter:
        %   note: this must match the algorithm
        fft_size = 2^nextpow2(N/4);
        
        % filter resampling mode:
        i_mode = 'pchip'; 

        % sampling rate:
        fs = 1;
        
        % DFT bin frequency step:
        bin_step = 2*fs/N;
        
        f0_rat_max = fs*p.f0_rat_max;
        f0_rat_min = p.f0_rat_min_bin*bin_step;
        if f0_rat_max == 0
            f0_rat_max = f0_rat_min;
        elseif f0_rat_max < f0_rat_min
            f0_rat_max = 1.05*f0_rat_min;
        end
        
        % select harmonic frequency:
        if isfield(p,'f0_rat') && p.f0_rat
            % fixed frequency mode:            
            rel_f = p.f0_rat*(1 + p.f0_rat_rnd*(2*rand-1));
            rel_f = max(rel_f,0);
            rel_f = min(rel_f,f0_rat_max);            
            f0_rat_max = 1;            
        else
            % randomized frequency mode:
            rel_f = rand(1);
        end
                
        % harmonic frequency:
        f0 = (f0_rat_max - f0_rat_min)*rel_f + f0_rat_min;
        
        %f0/fs
        
        % harmonic amplitude:
        a0 = p.f0_amp;
        
        % harmonic phase angle:
        ph = p.f0_phi + p.f0_phi_rnd*(2*rand(1) - 1);
                
        % 2*pi*time vector:
        tw = reshape([0:N-1]/fs*2*pi,[N 1]);
             
        
        % synthesize fundamental waveform:
        u = a0.*sin(tw*f0 + ph);
        
        % add some noise:
        u = u + randn(N,1)*p.rms_noise;
        
        % add random offset at the bitres level:
        u = u + (2*rand(1)-1).*a0.*2^-(p.bits-1);
        
        % round to ADC resolution:
        u = round(u/a0*2^(p.bits-1))/2^(p.bits-1)*a0;
               
        % generate filter freq. axis:
        ff = linspace(0,fs/2,fft_size+1)(:);
        
        % generate gain:
        fgain = ones(size(ff));
        
        % generate random phase gradient starting from 0:
        %  ####fixed complex results of negative slopes
        fphi_pwr = rand(1)*(2.5 - 0.5) + 0.5;
        fphi_max = p.ff_max_phi*(2*rand(1) - 1);
        fphi = sign(fphi_max)*linspace(0,abs(fphi_max),numel(ff)).^fphi_pwr/(abs(fphi_max)^(fphi_pwr-1));
        if any(isnan(fphi))
            fphi = zeros(size(fphi));
        end
        
        % generate random amplitude gradient starting from 0:
        %  ####fixed complex results of negative slopes
        famp_pwr = rand(1)*(2.5 - 0.5) + 0.5;
        famp_max = p.ff_max_amp*(2*rand(1) - 1);
        fgain = 1 + sign(famp_max)*linspace(0,abs(famp_max),numel(ff)).^famp_pwr/(abs(famp_max)^(famp_pwr - 1));
        if any(isnan(fgain))
            fgain = ones(size(fgain));
        end
        
        
        % make some window for analysis:
        w = flattop(N,8)(:);
        %  get window scaling factor:
        w_gain = mean(w);
        %  get window rms:
        w_rms = mean(w.^2).^0.5;
        
        % estiamte DC level:
        dc = mean(u.*w)/w_gain;
        
        % remove DC offset:
        u = u - dc;        
                        
        
        % apply filter:
        [uf,a,b, fff,ffg,ffp] = td_fft_filter(u, fs, fft_size, ff,fgain,fphi, i_mode);
        
        % make filtered and original waveform have identical length:
        u = u(a:b);
        
        % new samples count:
        M = numel(u);
        
        % make some window for analysis:
        w = flattop(M,8)(:);
        %  get window scaling factor:
        w_gain = mean(w);
        %  get window rms:
        w_rms = mean(w.^2).^0.5;
            
        
        % apply window to both signals:
        uc = [u uf].*w;
        
        
        % do FFT:
        U = fft(uc)(1:round(M/2),:)/M*2/w_gain;
        fh = [0:size(U,1)-1]'/M*fs;
        
        % fundamental component DFT bin:
        fid = round(f0/fs*M + 1);
        
        % filtered:     
        Uf = U(:,2);
        
        % reference:
        Ur = U(:,1);
                
        % phase correction value:
        fphib = interp1(fff,ffp,fh,'pchip','extrap');
        
        % amplitude correction value:
        fampb = interp1(fff,ffg,fh,'pchip','extrap');
        
        % apply filter correction to the reference signal:
        Ur = Ur.*fampb.*exp(j*fphib);       
               
                        
        
        % window size:
        w_size = 11;
        
        % not processed DFT bins:
        msk = [p.f0_rat_min_bin:numel(fh)];
                
        % DFT bins occupied by the harmonic
        h_bins = max((fid - w_size),1):min(fid + w_size,M);
        
        % remove harmonic bins from remaining list:
        msk = setdiff(msk,h_bins);
        msk = msk(msk <= N & msk > 1);
        
        % estimate noise levels for the removed harmonics components:
        Urns = interp1(fh(msk),abs(Uf(msk)),fh,'nearest','extrap');
        Ufns = interp1(fh(msk),abs(Ur(msk)),fh,'nearest','extrap');
        
        % estimate RMS noise from windowed spectrum:
        Ur_noise = sum(0.5*abs(Urns(w_size:end)).^2)^0.5/w_rms*w_gain;
        Uf_noise = sum(0.5*abs(Ufns(w_size:end)).^2)^0.5/w_rms*w_gain;
        
        
        % estimate RMS noise from windowed spectrum:
        Ur_rms = sum(0.5*abs(Uf(w_size:end)).^2)^0.5/w_rms*w_gain;
        Uf_rms = sum(0.5*abs(Ur(w_size:end)).^2)^0.5/w_rms*w_gain;
        
        
%         loglog(fh,abs(Uf))
%         hold on;
%         loglog(fh,abs(Ur),'r')
%         hold off;
        
        % amplitude deviation:
        dA(k) = abs(Ur(fid))/abs(Uf(fid)) - 1;
        
        % phase deviation:
        dP(k) = mod(arg(Ur(fid)) - arg(Uf(fid)) + pi, 2*pi) - pi;
        
        % noise gain:
        dN(k) = (Uf_noise/Ur_noise) - 1;              
        
        % rms gain:
        dR(k) = Uf_rms/Ur_rms - 1;
          
        
    endfor
    
    % detect worst cases:
    r.wA = single(max(abs(dA)));
    r.wP = single(max(abs(dP)));
    r.wN = single(max(abs(dN)));
    r.wR = single(max(abs(dR)));
    
    % estimate 95% uncertainty:
    if numel(dA)
        r.uA = single(est_scovint(dA,0));
        r.uP = single(est_scovint(dP,0));
        r.uN = single(est_scovint(dN,0));
        r.uR = single(est_scovint(dR,0));
    end
    
endfunction


