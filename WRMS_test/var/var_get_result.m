%% Get single result from result buffer. If no result ID specified, first one is taken.

function [r] = var_get_result(vr,res,id)

  if(~exist('id','var'))
    id = 1;
  end
  
  if(id>vr.res_n)
    error(['Result ID = ' int2str(id) ' out of range (measured count = ' int2str(vr.res_n) ')!']);
  end

  % return single result
  r = res{id};

endfunction