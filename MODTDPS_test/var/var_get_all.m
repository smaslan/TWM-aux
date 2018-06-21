function [outp] = var_get_all(par,vr,step,verbose)
% Original slow version of the 'var_get_all_fast()' function.
% Behaviour is the same, so see help of 'var_get_all_fast()'.
%
% License:
% --------
% This is part of VAR library for automatic multidim. variation of simulation parameters.
% (c) 2018, Stanislav Maslan, s.maslan@seznam.cz
% The script is distributed under MIT license, https://opensource.org/licenses/MIT 

  if(verbose)
    printf('Generating parameter combinations ... \r');
  end
  
  % get first combination
  [p,vr] = var_get_next(par,vr);
    
  % allocate full buffer
  outp = repmat(p,vr.var_n,1);
 
  if(vr.var_n>1)
    % generate rest of the combinations
    k = 2;
    da_end = 0;
    while(~da_end)
      % generate chunk
      pl = min(step,vr.var_n-k+1);
      [outp(k:k+pl-1),vr,da_end] = var_get_next(par,vr,step);
      k += pl;
      
      if(verbose)
        printf('Generating parameter combinations ... %3.0f%%  \r',100*outp(k-1)._vpid_/outp(k-1)._vpcnt_);
      end       
    end
          
  end
  
  if(verbose)
    printf('\n');
  end
    
  % convert to cells
  outp = num2cell(outp);

end