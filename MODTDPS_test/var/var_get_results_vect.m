%% Returns vector of results for defined parameter 'name'.
%% Returns results vector for first found vector parameter if 'name' undefined.
%% If multi-dimensional simulation has been used (multiple input vector parameters)
%% specify fixed indexes for rest of the vector parameters (other than 'name'):
%%  format: (..., 'param_1_name', param_1_index, 'param_2_name', param_2_index, ...)
%% 
%% Returns:
%%  r        - results cell vector
%%  v        - results cell vector coverted to struct of paramters (vectors or arrays)
%%  par_val  - vector of values of the variable parameter
%%  par_name - name of the variable parameter
%%  info_str - information about fixed parameters setup: 'param_1_name = param_1_value, param_2_name = param_2_value, ...'   

function [r,v,par_val,par_name,info_str] = var_get_results_vect(res,par,vr,varargin)

  if(~length(varargin))
    %% parameter name undefined, search first vector parameter
    pid = find(vr.par_n>1,1,'first');
    
    if(~isscalar(pid))
      % take first parameter if non of the parameters is vector
      pid = 1;
    end

  else
    %% user defined paramter name, check whether exist in paramters list
    pid = find(strcmp(vr.names,varargin{1}),1,'first');
    
    if(~isscalar(pid))
      % not found
      error('Invalid parameter name!');      
    end
              
  end
  
  %% no filtering yet
  fil_mult = zeros(1,vr.n);
  fil_ofs = zeros(1,vr.n);
  
  %% process find filters
  if(length(varargin)>2)
    
    %% function options: 'par_name_1', par_1_index, 'par_name_2', par_2_index, ... 
    varargin = {varargin{2:end}};
    
    % total optional paramters count
    n = length(varargin);
       
    % for each option
    for k = 1:2:n
      
      if(k<n)
        % find parameter in parameters list
        fid = find(strcmp(vr.names,varargin{k}),1,'first');
        if(~numel(fid))
          % parameter not found
          error(['Parameter "' inputname(1) '.' varargin{k} '" not found!']);          
        end
        
        if(varargin{k+1}<1 || varargin{k+1}>vr.par_n(fid))
          % parameter index out of range
          error(['Parameter "' inputname(1) '.' varargin{k} '" index out of range (' int2str(varargin{k+1}) ' of ' int2str(vr.par_n(fid)) ')!']);          
        end
        
        % filter multiplier
        fil_mult(fid) = prod([1;vr.par_n(1:fid-1)]);
        
        % filter offset
        fil_ofs(fid) = varargin{k+1}-1;
      end      
    end
  end
  
  % finish filter
  ofs = fil_mult*fil_ofs';
  
  %% get results diplacement in results buffer for selected vector of parameters
  idx_step = prod([1;vr.par_n(1:pid-1)]);
  
  %% vector length
  idx_rng = vr.par_n(pid);
  
  %% return results vector
  r = {res{(0:idx_rng-1)*idx_step+1 + ofs}};
  
  %% convert vector of result structs to struct of result vectors/arrays (hardcore implementation)
  v = vectorize_structs_elements(r);
  
  %% return variable parameter values
  par_val = getfield(par,vr.names{pid});
  
  %% return name of the variable parameter
  par_name = vr.names{pid};
     
  %% build info string for the fixed parameters
  % for each parameter
  info_str = '';
  for k = 1:vr.n
    if(vr.par_n(k)>1 && k~=pid)
      
      % only for vector parameters and not the variable parameter (var_name)
      if(numel(info_str))
        info_str = [info_str ', '];
      end
      
      % get parameter value
      val = getfield(par,vr.names{k})(fil_ofs(k)+1);
      
      % add to info string
      info_str = [info_str vr.names{k} '(' int2str(fil_ofs(k)+1) ') = ' num2str(val)];
            
    end
  end
   
endfunction
