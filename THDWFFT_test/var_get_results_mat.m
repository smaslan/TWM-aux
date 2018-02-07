%% Returns 2D martixes of results from multidimensional space of results.
%% Takes only scalars elements from results structure 'res'. Combines them 
%% into 2D matrixes. The matrixes are stores into output structure 'v'. 
%% Vector and matrix elements of 'res' structure are ignored.
%%
%% Parameters:
%%  res      - results cell vector
%%  par      - struct of simulation parameters (scalars and vectors)
%%  vr       - initialized control structure (see var_init())
%%  varargin - list of 'par' struct vector names:
%%    1)   'x_vector_name' - X vector name
%%    2)   'y_vector_name' - Y vector name
%%    3-?) 'vector_name_1',vector_1_index,'vector_name_2',vector_2_index
%% 
%% Returns:
%%  v        - struct of 2D matrices of results
%%  x_val    - values of: 'par.x_vector_name'
%%  y_val    - values of: 'par.y_vector_name'
%%  x_name   - name of X vector
%%  y_name   - name of Y vector
%%  info_str - string with list of the rest of vector elements of 'par' structure
%%           - 'elem_1(index_1) = value_at_index_1, elem_2(index_2) = value_at_index_2, ...'

function [v,x_val,y_val,x_name,y_name,info_str] = var_get_results_mat(res,par,vr,varargin)
  
  if(length(varargin)<2)
    error('Not enough parameters!');
  endif
    
  % get X vector length
  [len,x_val] = var_get_par_by_name(vr,par,varargin{1});
  x_name = varargin{1};
   
  % build set of paramter filters for Y vector
  varargin = {varargin{2} varargin{1} 1 varargin{3:end}};
         
  %%%% for every X vector item %%%%
  for k = 1:len
    %% build set of paramter filters for Y vector
    varargin{3} = k;
        
    %% get Y vector
    [r,v,y_val,y_name,info_str] = var_get_results_vect(res,par,vr,varargin{:});
       
    %% allocate buffers for every struct item
    if(k == 1)            
      % struct item names
      names = fieldnames(v);
      
      % struct items count
      pn = length(names);
                  
      % allocate object for every struct field
      mat = cell(pn,1);
      matn = cell(pn,1);
      n = 0;
      for m = 1:pn
        
        % get field
        fld = getfield(v,names{m});
        
        if(isvector(fld))
          % only for vectors
          mat{++n} = zeros(length(fld),len);
          matn{n} = names{m};
        end      
      endfor
      mat = {mat{1:n}};
      matn = {matn{1:n}};      
      pn = n; 
    endif
        
    %% combine vectors into matrix for every struct element
    for n = 1:pn
      % get struct files 
      fld = getfield(v,matn{n});
      
      % add to the matrix
      mat{n}(:,k) = fld;
    endfor
           
  endfor
  
  %% combine 2D matrices into structure
  v = struct();
  for n = 1:pn
    v = setfield(v,matn{n},mat{n});  
  endfor
        
  %% remove fist Y parameter name from info_str (hardcore)
  strrem = [varargin{2} '(' int2str(len) ') = '];
  cid = findstr(info_str,strrem);
  if(numel(cid))
    eid = cid + find(info_str(cid:end)==',',1,'first') - 1;
    if(numel(eid))
      if(cid==1)
        info_str = info_str(eid+2:end);
      else
        info_str = [info_str(1:cid-3) info_str(eid:end)];
      endif
    else
      info_str = info_str(1:cid-3);
    endif  
  endif
    
endfunction