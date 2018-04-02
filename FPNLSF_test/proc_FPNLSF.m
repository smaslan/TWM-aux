function [r] = proc_FPNLSF(p)
    
    
    % do monte-carlo
    dpx = [];
    dox = [];
    dAx = [];
    dfx = [];
    for r = 1:p.R
    
        % sampling rate [Hz]:
        fs = p.f0*(p.fs_rat + rand(1) - 0.5);
        
        % samples count:
        N = round(fs/p.f0*(p.f0_per + rand(1) - 0.5));
        
        % haropnic components:
        fh = [];
        fh(1,:) = p.f0:p.f0:floor(0.4*fs);        
        if p.sfdr_nc
            fh(2:end) = fh(2:end) + (rand(size(fh(2:end))) - 0.5)*p.f0;
        end
        
        % randomize harmonic spurr values:        
        A = [1 rand(size(fh(2:end)))*10^(-p.sfdr/20)];            
                
        % generate time vector:
        t = [];
        t(:,1) = [0:N-1]/fs + p.jitt*randn(1,N);
        
        % generate random phases of harmonics:
        phi = rand(size(fh))*2*pi;
        ph0 = phi(1);
        
        % synthesize waveform:
        u = sum(A.*sin(2*pi*t.*fh + ph0),2);
        
        % add some random offset: 
        ofs = rand(1)/2^(p.bits-1);
        u = u + ofs;
                
        
        for tr = 1:p.cfg.max_try
            % randomize initial guess (from second try):
            rand_p = randn(1)*0.01*pi*(tr > 1);
            rand_f = p.cfg.max_f_dev*randn(1)*fh(1)*(tr > 1);
            rand_o = 0.001*randn(1)*(tr > 1);
            
            % randomize offeset:
            ux = u + rand_o;
            
            % round to ADC resolution:
            ux = round(ux*2^(p.bits-1))/2^(p.bits-1);
            
            % estimate initial phase:
            phi_zc = phase_zero_cross(ux) + rand_p;
            if isinf(phi_zc)
                phi_zc = randn(1)*pi;
            end                               
            
            % fit waveform:            
            [Ax, fx, phx, ox] = FPNLSF(t,ux,fh(1)*(1+rand_f),0, phi_zc);
            ox = ox - rand_o;
                        
            if ~isinf(fx) && abs(fx/fh(1)-1) < p.cfg.max_f_dev
                % result possibly ok - leave
                fail = 0; 
                break;
            elseif tr == p.cfg.max_try
                disp('Warning: No convergence even after all retrires! Dunno what to do now...');
                fail = 1;
            end
            % retry because we got no convergence or too high phase deviation
                     
        end
        
        
        if ~fail
            % calculate phase deviation (wrapped to +-pi):
            dphx = (phx - ph0);
            dpx(end+1) = mod(dphx + pi,2*pi) - pi;
                                    
            % deviation of ampl.:
            dAx(end+1) = Ax - 1;
            
            % deviation of freq.:
            dfx(end+1) = fx/fh(1) - 1;
            
            % deviation of offset.:
            dox(end+1) = ox - ofs;
            
            %fprintf('o = %8g  p = %8g\n',dox(end),dpx(end));
        end
    
    end
    

    % 
    ids = ~(scores(dfx,2,0.9) | scores(dAx,2,0.9) | scores(dpx,2,0.9) | scores(dox,2,0.9));
    
    m_dfx = mean(dfx(ids));
    s_dfx = std(dfx(ids));    
    m_dAx = mean(dAx(ids));
    s_dAx = std(dAx(ids));    
    m_dpx = mean(dpx(ids));
    s_dpx = std(dpx(ids));
    m_dox = mean(dox(ids));
    s_dox = std(dox(ids));
    
    
    r = struct();
    % info
    r.num = sum(ids);
    % return mc lists:
    r.s_dpx = s_dpx;
    r.s_dAx = s_dAx;
    r.s_dfx = s_dfx;
    r.s_dox = s_dox;
    r.m_dpx = m_dpx;
    r.m_dAx = m_dAx;
    r.m_dfx = m_dfx;
    r.m_dox = m_dox;
    

end


