%% Get next combination of parameters. Set 'count' to 0 to return all combinations.

function [outp,vr,da_end] = var_get_next(par,vr,count)   
  
  % return single combinarion if count undefined
  if(nargin()<3)
    count = 1;
  endif
  % return all combinations for zero 'count'
  if(count==0)
    count = vr.var_n;
  endif
     
  %% prepare list of variable parameters
  % vector parameter id's
  vids = find(vr.par_n>1); 
  
  % vector paramter names
  v_names = {vr.names{vids}};
  
  % vector paramter indexes
  v_cnt = vr.par_cnt(vids);
  
  % vector paramter lengths
  v_max = vr.par_n(vids);
  
  % vector parametrs count
  v_n = length(vids);
  
  % cell array of vetor paramters values
  v_vals = cellfun(@getfield,repmat({par},1,v_n),v_names,'UniformOutput',false');
      
  % generate output prototype
  outprot = par;
  for p = 1:v_n
    outprot = setfield(outprot,v_names{p},0);
  endfor
  
  % allocate output vector
  outp = repmat(outprot,min(count,vr.var_n),1);
  
  %% return up to 'count' combinationations
  da_end = 0;
  k = 0;
  while(~da_end && k<count)
      
    % copy current set of parameters
    outp(++k) = outprot;
    for p = 1:v_n
      outp(k) = setfield(outp(k),v_names{p},v_vals{p}(v_cnt(p)));
    endfor
    % store combination id
    outp(k)._vpid_ = ++vr.var_id;
             
    % prepare next combination of parameters
    for p = 1:v_n
      if(v_cnt(p)>=v_max(p))
        v_cnt(p) = 1;
        da_end = (p==v_n);
      else
        v_cnt(p)++;
        break;
      endif          
    endfor
  endwhile
  
  % store work variables back to the 'vr' struct
  vr.par_cnt(vids) = v_cnt;
        
  % reduce allocated cells vector to actual size
  if(k==1)
    outp = outp(1);
  else
    outp = outp(1:k);
  endif
    

endfunction
