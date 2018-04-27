function [r] = proc_WRMS(p)

    warning('off');

    dP = [];
    dS = [];
    dQ = [];
    dI = [];
    dIn = [];
    for k = 1:p.cycles

        % normalized fundamental frequency:
        f0 = 1;
    
        % samples per period +- 1:
        fs = p.fs_rat - (2*rand(1) + 1);
        
        % channel A amplitude:
        A = 1;
        
        % channel B amplitude +- some gain noise:
        B = A*p.ab_rat.*(1 + p.ab_rat_rnd*(2*rand(1) - 1));
        
        % interchannel phase shift +- random phase-step range: 
        dph = p.ab_phi + p.ab_phi_rnd*(2*rand(1) - 1);    
        
        % generate absolute phases (1,channel):
        ph = rand(1)*2*pi + [0 dph];
        
        % generate samples count +- some length variation:
        N = round(p.f0_per*(1 + 0.10*(2*rand(1) - 1))*fs);
                
        % 2*pi*time vector:
        tw = reshape([0:N-1]/fs*2*pi,[N 1]);
                
        % synthesize fundamental waveforms:
        u_clean = [A B].*sin(tw + ph);
        u = u_clean;
        
        % pregenerate noise:
        rnd_noise = randn(N,2)*p.rms_noise;
        
        % pregenerate random offset:
        rnd_dc = randn(1,2).*[A B].*[2^-(p.a_bits-1) 2^-(p.b_bits-1)];
        
        % generate window:
        w = blackmanharris(N,'periodic');
        
        % calculate inverse RMS of the window (needed for scaling of the result): 
        W = mean(w.^2)^-0.5; 
                        
        
        for s = 1:p.spurr_mode+1
        
            if s > 1
                % -- add spurr in spurr mode and second pass only:
                
                % restore clean harmonics:
                u = u_clean;
            
                % spurr frequency: 
                f_spurr = p.s_freq*(fs/2 - f0) + f0;
                            
                % add spurr waveform to B channel:        
                u(:,2) = u(:,2) + p.s_amp*sin(tw*f_spurr + rand(1)*2*pi);
                
                % spurr I rms:
                Is = p.s_amp*2^-0.5;              
            else
                % spurr I rms:
                Is = 0;
            end
            
            % add some noise:
            u = u + rnd_noise;
            
            % add random offset:
            u = u + rnd_dc;
            
            % round to ADC resolution:
            u(:,1) = round(u(:,1)*2^(p.a_bits-1))/2^(p.a_bits-1);
            u(:,2) = round(u(:,2)/B*2^(p.b_bits-1))/2^(p.b_bits-1)*B;
                                  
            % apply window to channels:
            uw = u.*w;
            
            % estimate DC offset:
            dc = mean(uw,1).*W^2;
            
            % remove DC offset:
            u = u - dc;
            
            % apply window to channels:
            uw = u.*w;
            
            % eval reference values in freq. domain:
            %  with no spurr:
            Pr = (0.5*A*B*cos(dph));
            Qr = abs(0.5*A*B*sin(dph));
            Sr = 0.5*A*B;
            Ir = (2^-0.5*B);
                        
            % eval U/I rms values in time domain:
            UI = W.*mean(uw.^2,1).^0.5;
            Ui = UI(1); 
            Ii = UI(2);
            
            % remove spurr from rms current from time domain:
            %  note: we are interested to error of the fundamental caused by the spurr
            Ii = (Ii^2 - Is^2).^0.5;
            
            % apparent power in time domain:
            Si = Ui*Ii;
            
            % eval power in time domain:
            Pi = W^2*mean(prod(uw,2),1);
            
            % reative power in time domain:
            Qi = (Si.^2 - Pi.^2).^0.5;
            
            if s == 1
                % -- single tone mode
            
                % power deviation:
                dP(k) = abs(Pi - Pr)/Sr;            
                % apparent power deviation:
                dS(k) = abs(Si - Sr)/Sr;            
                % reactove power deviation:
                dQ(k) = abs(Qi - Qr)/Sr;            
                % rms current deviation:
                dI(k) = abs(Ii - Ir)/Ir;
                
                % store spurr free values for next pass:
                Pi0 = Pi;
                Si0 = Si;
                Qi0 = Qi;
                Ii0 = Ii;
                
            else
                % -- spurr mode
                %  note: store devaitions spurr - no-spurr
                
                % power deviation:
                dP(k) = abs(Pi - Pi0)/Sr;            
                % apparent power deviation:
                dS(k) = abs(Si - Si0)/Sr;            
                % reactove power deviation:
                dQ(k) = abs(Qi - Qi0)/Sr;            
                % rms current deviation:
                dI(k) = abs(Ii - Ii0)/Ir;

            end
        
        end
        
    endfor
    
    % detect worst cases:
    r.dP = single(max(dP));
    r.dS = single(max(dS));
    r.dQ = single(max(dQ));
    r.dI = single(max(dI));
    
endfunction


