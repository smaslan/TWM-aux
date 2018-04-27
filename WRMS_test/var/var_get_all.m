% Get all paramter combinations.

function [outp] = var_get_all(par,vr,step,verbose)

  if(verbose)
    printf('Generating parameter combinations ... \r');
  endif
  
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
      endif          
    endwhile
          
  endif
  
  if(verbose)
    printf('\n');
  endif
    
  % convert to cells
  outp = num2cell(outp);

endfunction