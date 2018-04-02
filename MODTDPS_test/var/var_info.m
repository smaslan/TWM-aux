%% Prints simulation progress.

function [] = var_info(par,verbose)

  %% simulation info
  if(verbose)
    disp(['=== SIMULATION: cycle ' int2str(par._vpid_) ' of ' int2str(par._vpcnt_) ' ===']);
  end
    
  if(verbose>1)
    %% display current parameters info
    
    % get struct member names
    names = fieldnames(par);
    
    % find maximum parameter name length
    name_len = max(cellfun(@length,names));
            
    % print value for each parameter 
    for p = 1:length(names)
      str = repmat(' ',1,name_len);
      str(1:length(names{p})) = names{p};
      
      val = getfield(par,names{p});
      if(isscalar(val) && ~isstruct(val) && ~strcmp(names{p},'_vpid_') && ~strcmp(names{p},'_vpcnt_'))
        disp(['  ' str ' = ' num2str(val)]);
      end          
    end
    
    disp('');
  end
  
  

endfunction