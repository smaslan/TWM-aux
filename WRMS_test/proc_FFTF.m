function [r] = proc_WRMS(p)

    warning('off');

    
    dA = [];
    dP = [];
    dN = [];
    dR = [];
    for k = 1:p.cycles
    
        % random size from 2^(p.N_pow-1)+1 to 2^(p.N_pow) 
        N = 2^p.N_pow - round((2^(p.N_pow-1) - 1)*rand(1));
        
        % size of the FFT filter:
        %   note: this must match the algorithm
        fft_size = 2^nextpow2(N/4);
        
        % filter resampling mode:
        %  note: we will generate the filter at its expected size, so we can use 'nearest'
        i_mode = 'nearest'; 

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
                
        % harmonic frequency:
        f0 = (f0_rat_max - f0_rat_min)*rand(1) + f0_rat_min;
        
        %f0/fs
        
        % harmonic amplitude:
        a0 = p.f0_amp;
        
        % harmonic phase angle:
        ph = p.f0_phi_rnd*(2*rand(1) - 1);
                
        % 2*pi*time vector:
        tw = reshape([0:N-1]/fs*2*pi,[N 1]);
             
        
        % synthesize fundamental waveform:
        u = a0.*sin(tw*f0 + ph);
        
        % add some noise:
        u = u + randn(N,1)*p.rms_noise;
        
        % add random offset at the bitres level:
        u = u + randn(1).*a0.*2^-(p.bits-1);
        
        % round to ADC resolution:
        u = round(u/a0*2^(p.bits-1))/2^(p.bits-1)*a0;
               
        % generate filter freq. axis:
        ff = linspace(0,fs/2,fft_size+1)(:);
        
        % generate gain:
        fgain = ones(size(ff));
        
        % generate random phase gradient starting from 0:
        fphi_pwr = rand(1)*(3 - 0.1) + 0.1;
        fphi_max = p.ff_max_phi*(2*rand(1) - 1);
        fphi = linspace(0,fphi_max,numel(ff)).^fphi_pwr/(fphi_max^(fphi_pwr-1));
        if any(isnan(fphi))
            fphi = zeros(size(fphi));
        end
        
        % generate random amplitude gradient starting from 0:
        famp_pwr = rand(1)*(3 - 0.1) + 0.1;
        famp_max = p.ff_max_amp*(2*rand(1) - 1);
        fgain = 1 + linspace(0,famp_max,numel(ff)).^famp_pwr/(famp_max^(famp_pwr - 1));
        if any(isnan(fgain))
            fgain = ones(size(fgain));
        end
                        
        
        % apply filter:        
        [uf,a,b] = td_fft_filter(u, fs, fft_size, ff,fgain,fphi, i_mode);
        
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
        U = fft(uc)(1:round(M/2),:);
        fh = [0:size(U,1)-1]'/M*fs;
        
        % fundamental component DFT bin:
        fid = round(f0/fs*M + 1);
        
        % reference vector:
        Ur = U(fid,1);
        
        % mathing frequency of the DFT bin:
        f0b = (fid - 1)/M*fs;
        
        % phase correction value:
        fphib = interp1(ff,fphi,f0b,'pchip','extrap');
        
        % amplitude correction value:
        fampb = interp1(ff,fgain,f0b,'pchip','extrap');
        
        % apply filter correction to the reference signal:
        Ur = Ur.*fampb.*exp(j*fphib);       
        
        % filtered vector:     
        Uf = U(fid,2);
        
        
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
        Urns = interp1(fh(msk),abs(U(msk,1)),fh,'nearest','extrap');
        Ufns = interp1(fh(msk),abs(U(msk,2)),fh,'nearest','extrap');
        
        % estimate RMS noise from windowed spectrum:
        Ur_noise = sum(0.5*Urns.^2)^0.5/w_rms*w_gain;
        Uf_noise = sum(0.5*Ufns.^2)^0.5/w_rms*w_gain;
        
        % estimate RMS noise from windowed spectrum:
        Ur_rms = sum(0.5*abs(U(3:end,1)).^2)^0.5/w_rms*w_gain;
        Uf_rms = sum(0.5*abs(U(3:end,2)).^2)^0.5/w_rms*w_gain;
        
        
        %loglog(fh,abs(U(:,1)))
        %hold on;
        %loglog(fh,abs(U(:,2)))
        %hold off;
               
        
        % amplitude deviation:
        dA(k) = abs(Ur)/abs(Uf) - 1;
        
        % phase deviation:
        dP(k) = mod(arg(Ur) - arg(Uf) + pi, 2*pi) - pi;
        
        % noise gain:
        dN(k) = (Uf_noise/Ur_noise) - 1;              
        
        % rms gain:
        dR(k) = Uf_rms/Ur_rms - 1;
          
        
    endfor
    
    % detect worst cases:
    r.dA = single(max(abs(dA)));
    r.dP = single(max(abs(dP)));
    r.dN = single(max(abs(dN)));
    r.dR = single(max(abs(dR)));
    
endfunction


