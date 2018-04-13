% Get all parameter combinations - faster version.

function [outp] = var_get_all_fast(par,vr,step,verbose)

    if(verbose)
        printf('Generating parameter combinations ... \r');
    endif
    
    % total combinations:
    N = prod(vr.par_n);
        
    % create prototype:
    p_prot = par;
    
    c_prot = {}; 
    % replace all vector parameters by default scalars:
    for v = 1:numel(vr.names)                
        if vr.par_n(v) > 1
            p_prot = setfield(p_prot,vr.names{v},0); 
        end
        c_prot(v,1) = getfield(p_prot,vr.names{v});
    end
    
    vidid = find(strcmpi(vr.names,'_vpid_'));
    
    % vector variables count:
    vn = sum(vr.par_n > 1);
    
    % vector variables ids:
    vids = [1:numel(vr.par_n)](vr.par_n > 1);

    % load vectors:
    vars = {};        
    for v = 1:vn
    
        % get variable vector:
        v_val = getfield(par,vr.names{vids(v)});
        
        % reshape the vector to v-dim:
        dimn = [(eye(vn)(v,:))*(numel(v_val)-1) + 1];
        vars{v} = reshape(v_val,dimn);
    
    end
    
    % create combinations matrix (combination,variable):
    vc_lists = [];
    for v = 1:vn    
        v_prod = 1;
        for d = 1:vn
            if d == v
                v_part = vars{d};                
            else
                v_part = ones(size(vars{d}));
            end
            v_prod = bsxfun(@times,v_part,v_prod);
        end        
        vc_lists(:,v) = v_prod(:);
    end
    
    
    % initialize vector of parameters:
    c_list = repmat({c_prot},[1 N]);
    
    % store variable combinations:
    tot = 0;
    while tot < N             
        todo = min(N - tot,step);
        for k = (tot+1):(tot+todo)
            for v = 1:vn
                c_prot{vids(v)} = vc_lists(k,v);
                
            end
            c_prot{vidid} = k;
            c_list{k} = c_prot;
        end
        tot += todo;
        if(verbose)
            printf('Generating parameter combinations ... %3.0f%%  \r',100*tot/N);
        endif   
    end
  
    if(verbose)
        printf('\n');
    endif
    
    % make cell array of parameter structs:
    outp = cellfun(@cell2struct,c_list,{vr.names},{1});
    
    % convert to cells
    outp = num2cell(outp);
    
endfunction