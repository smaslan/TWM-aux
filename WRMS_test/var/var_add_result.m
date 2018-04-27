%% Adds new result 'newres' to results buffer 'res'.

function [res,vr] = var_add_result(res,vr,newres)
  
  % add result to a buffer
  res{++vr.res_n} = newres;
    
end