function [r] = proc_MODTDPS(p)

    warning('off');

    if p.is_rect
        wshape = 'rect';
    else
        wshape = 'sine';
    end
    
    % do monte-carlo
    dofs = [];
    df0  = [];
    dfm  = [];
    dA0  = [];
    dmod = [];
    for r = 1:p.R_max
    
        % sampling rate [Hz]:
        fs = p.fsf0_rat + rand(1) - 0.5;
        
        % carrier frequency [Hz]:
        f0 = 1;
        
        % modulating signal freq [Hz]:
        fm = p.fmf0_rat*(1 + 0.02*(rand(1)-0.5));
        
        % samples count:
        N = round(fs/fm*(p.fm_per + rand(1) - 0.5));
        
        % carrier amplitude:        
        A0 = 1;            
        
        % modulating amplitude:        
        Am = min(max(A0*p.modd*(1 + 0.05*(rand(1)-0.5)),0.01),0.99);
        
        % pk-pk amplitude:
        App = A0*(1 + p.modd);
        
        % absolute ADC resolution:
        Abr = App*2^-p.bits;
        
        % randomize phases:
        ph0 = rand(1)*2*pi;
        phm = rand(1)*2*pi;
        
        % randomize offset:
        ofs = A0*p.ofs*randn(1);
                
        % generate time vector:
        t = [];
        t(:,1) = ([0:N-1]/fs + p.jitt*randn(1,N))*2*pi;
                
        % synthesize signal:
        if p.is_rect
            u = ofs + sin(t*f0 + ph0).*(A0 + Am*(0.5 - (mod(t*fm + phm,2*pi) > pi)));
        else
            u = ofs + sin(t*f0 + ph0).*(A0 + Am*sin(t*fm + phm));
        end
                
        
        % -- add THD/spurrs
        % haropnic components:
        fh = [];
        fh(1,:) = (2*f0):f0:floor(0.4*fs);        
        if p.sfdr_nc
            fh = fh + (rand(size(fh)) - 0.5)*f0;
        end
        % randomize harmonic spurr values:        
        Ah = rand(size(fh))*10^(-p.sfdr/20);
        % generate random phases of harmonics:
        phh = rand(size(fh))*2*pi;   
        % add harmonic distortion:
        u = u + sum(Ah.*sin(t.*fh + phh),2);
        
        % round to ADC resolution: 
        u = round(u/Abr)*Abr;
        
        % evaluate algortithm:
        try 
            [me, dcx,f0x,A0x, fmx,Amx,phmx] = mod_tdps(fs, u, wshape, p.comp_err);
                                   
            modx = Amx/A0x;            
            modd = Am/A0;
            
            if ~isnan(f0x) && ~isnan(fmx) && ~isnan(modx) && ~isnan(modd) && ~isinf(f0x) && ~isinf(fmx) && ~isinf(modx) && ~isinf(modd)            
                dofs(end+1) = dcx - ofs;
                df0(end+1)  = f0x/f0 - 1;
                dfm(end+1)  = fmx/fm - 1;
                dA0(end+1)  = A0x/A0 - 1;
                dmod(end+1) = modx - modd;
            end
        end
        
        if numel(dmod) > p.R
            break;
        end
        
    end
        
    r = struct();
    % return mc lists:
    r.l_dofs = single(dofs);
    r.l_df0 = single(df0);
    r.l_dfm = single(dfm);
    r.l_dA0 = single(dA0);
    r.l_dmod = single(dmod);
    r.num = numel(dofs);
    
    if r.num > 2
        r.s_dofs = single(est_scovint(dofs,0));
        r.s_df0 = single(est_scovint(df0,0));
        r.s_dfm = single(est_scovint(dfm,0));
        r.s_dA0 = single(est_scovint(dA0,0));
        r.s_dmod = single(est_scovint(dmod,0));    
    end
    
    r.m_dofs = single(mean(dofs));
    r.m_df0 = single(mean(df0));
    r.m_dfm = single(mean(dfm));
    r.m_dA0 = single(mean(dA0));
    r.m_dmod = single(mean(dmod));

end


