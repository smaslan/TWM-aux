%% Initialize automatic parameter variation algorithm.

function [vr,par] = var_init(par)
  
  % add paramter combination ID into the parameters structure
  par._vpid_ = 0;
  
  % add total paramter combinations count into the parameters structure
  par._vpcnt_ = 0;
      
  % get input paramter names
  vr.names = fieldnames(par);
  
  % get parameters count
  vr.n = length(vr.names);
  
  % create variation counters for each paramter
  vr.par_cnt = ones(1,vr.n);
      
  % get parameter types and lengths (1 for scalar, N for vector)
  vr.par_n = cellfun(@length,cellfun(@getfield,repmat({par},length(vr.names),1),vr.names,'UniformOutput',false));
  
  % get total variations count
  vr.var_n = prod(vr.par_n);
  par._vpcnt_ = vr.var_n;
  
  % no paramter combinations generated yet
  vr.var_id = 0;
  
  % no results measured yet
  vr.res_n = 0;
  
endfunction

