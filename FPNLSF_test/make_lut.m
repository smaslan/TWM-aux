function [lut] = make_lut(res,p,vr,ax,qu)
% Simple generator of multidim lookup table (LUT) based on the data from 
% 'var_' functions (automatic variation of N-dim parameter space).
%
% The point of the VAR lib is to generate all possible combinations
% of input vector parameters, then calculate ANYTHING for every combination 
% and return the n-dim space of results, one dim for each input parameter.
% It was designed mainly for uncertainty and sensitivity analysis of 
% algorithms.  
%
% This function takes the results from the VAR lib and converts them
% to the compact LUT. Complement function 'interp_lut()' will interpolate
% in the LUT. It was designed for an uncertainty estimator based on the 
% precalculated data, but it is usable for anything. 
%
%
% Syntax:
%   unc = make_lut(res,p,vr,ax,qu)
%
% Inputs:
%   res - cell array of results from VAR lib.
%   p   - input parameters used for generation of 'res'.
%   vr  - control struct of the VAR lib.
%   ax  - struct of axes (parameters), one element per axis:  
%           ax.axis_1... - config of axis 'axis_1'
%           ax.axis_2... - config of axis 'axis_2'
%           ...
%         The 'axis_?' must match to the parameter names in the 'p'
%         and there must be one axis for each parameter 'p' of non-unity
%         size!
%         Each axis structure contains following:
%           axis.min_ovr - minimum allowed overrange, e.g. 0.95 means the
%                          interpolator function won't generate 
%                          out-of-range event if value >=95% of min axis value
%           axis.max_ovr - minimum allowed overrange, e.g. 1.05 means the
%                          interpolator function won't generate 
%                          out-of-range event if value >=105% of max axis value
%           axis.min_lim - min-out-of-range action ('error' to generate error, 
%                          or 'const' to limit interpolating value to min
%                          available value of axis)
%           axis.max_lim - max-out-of-range action ('error' to generate error, 
%                          or 'const' to limit interpolating value to max
%                          available value of axis)
%           axis.scale   - expected scaling of the axis. If the values 
%                          are distributed linearly, the interpolator
%                          will probably work ok with 'lin' option. If
%                          the data are displaced in log-scale, the option
%                          'log' may be better. The interpolator will apply
%                          log10() before interpolating.
%   qu - struct of quantities to be stored to the LUT:
%          val.quantity_1... - config of quantity 'quantity_1'
%          val.quantity_2... - config of quantity 'quantity_2'
%        Note the function will work only for scalar quantities!
%        Each quantity MAY contain following:
%          quantity.scale - scaling of the quantity before interpolation.
%                           'lin' or 'log' allowed. 'log' will apply 'log10(q)'
%                           before interpolation and then '10^q'. This may be useful
%                           when the quantity changes over many decades. Default
%                           is 'lin'.
%          quantity.mult - optional multiplier of the quantity after interpolation.
%                          default is 1.0.   
%
% Returns:
%   lut - lookup table structure:
%           lut.ax        - copy of axes 'ax'
%           lut.ax.values - copy of axis values from 'p'
%           lut.qu        - copy of quantity setups 'qu'
%           lut.qu.quantity_1.data - N-dim matrix of quantity values
%                                  - order of axes matches the order in 'ax'.
%                                  - format description see below
%           lut.qu.quantity_1.data_scale - scaling factor for 'data'         
%           lut.qu.quantity_1.data_offset - offset for 'data'
%           lut.qu.quantity_1.data_mode - format string, now only 'log10u16' 
%           ...
%
% 'log10u16' format: the function applies log10(q) to the quantity, then
% detects range qmin and qmax, then scales the quantity to 16bit unsinged integer
% and stores the result to the structure. To decode:
%   decode_data = 10.^(qu.data*qu.data_scale + qu.data_offset);
% Note it was designed for uncertainties where the accuracy is not critical.
% It of course cannot handle negative or zero values! 
%
% License:
% --------
% (c) 2018, Stanislav Maslan, smaslan@cmi.cz
% The script is distributed under MIT license, https://opensource.org/licenses/MIT

    % get required axes:
    ax_names = fieldnames(ax);
    A = numel(ax_names);
    
    % existing axes:
    v_names = vr.names(vr.par_n > 1);
    V = numel(v_names);
    
    % existing axes sizes:
    adims = vr.par_n(vr.par_n > 1);
    
    % init LUT:
    lut.ax = struct();
    
    % check is required axes exist in the data:
    m = A;
    for k = 1:V
        % look for the axis in the calculated data:
        vid = strcmpi(v_names{k},ax_names);
        if ~any(vid)
            error(sprintf('Axis ''%s'' not found in the axes configuration list ''ax''!',v_names{k}));            
        end
        m--;
        
        % get axis:
        arec = getfield(ax,v_names{k});
        
        % store axis values:
        arec = setfield(arec,'values',getfield(p,v_names{k}));
        
        % store axis record:
        lut.ax = setfield(lut.ax,v_names{k},arec);
                        
    end
    % check if all avaliable axes are assinged:
    if m
        error('Axis configuration ''ax'' contains at least one more axis that available in the calculated data!');
    end
       
    
    
    % get required quantity names:
    q_names = fieldnames(qu);
    Q = numel(q_names);
    
    % copy quantity setup to the LUT:
    lut.qu = qu; 
    
    % --- process each quantity:
    for k = 1:Q
    
        % quantity name:
        name = q_names{k};
        
        disp(sprintf('Processing quantity ''%s'' (%d of %d)',name,k,Q));
    
        % get the quantity element from each cell:
        try 
            data = cellfun(@getfield,res,{name});
        catch
            error(sprintf('Quantity ''%s'' not found in the calculated data!',name));
        end
               
        % shape the data to N-dim array matching the data orientation:
        data = reshape(data,adims);
        
        % -- compression:
        % convert to log-scale:
        c_data = log10(data);
        
        % detect range:
        qu_min = min(c_data(:));
        qu_max = max(c_data(:));
        
        % rescale to uint16:
        int_max = 2^16 - 1;                
        c_data = uint16(round(int_max*(c_data - qu_min)/(qu_max - qu_min)));
                
        % get quantity record:
        qu_rec = getfield(lut.qu,name);
        
        % store data and scaling factors:
        qu_rec.data = c_data;
        qu_rec.data_scale = (qu_max - qu_min)/int_max;
        qu_rec.data_offset = qu_min;
        qu_rec.data_mode = 'log10u16';
        
        % try to decode back to original values:
        b_data = 10.^(double(c_data)*qu_rec.data_scale + qu_rec.data_offset);
        % calculate deviation of the compressed data:
        dev = max(abs(b_data./data-1)(:));
        
        % store the quantity back to LUT:
        lut.qu = setfield(lut.qu,name,qu_rec);
            
    end         

end