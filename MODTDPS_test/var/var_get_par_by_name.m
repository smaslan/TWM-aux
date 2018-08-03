%% Get parameter length an values by its name

function [len,par_val] = var_get_par_by_name(vr,par,name)

  % find parameter id
  pid = find(strcmp(vr.names,name),1,'first');
  
  if(~numel(pid))
    % paramter not found
    error(['Paramter ''' name ''' is not member of paramters structure!']);  
  endif
  
  % return length of the parameter
  len = vr.par_n(pid);
  
  % return vector of parameters
  par_val = getfield(par,vr.names{pid});

end