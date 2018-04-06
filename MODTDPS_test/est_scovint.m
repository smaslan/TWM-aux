function [unc] = est_scovint(x,x0)    
    [nix,sql,sqr] = scovint(x,0.95,x0);
    if isempty(sql)
        [nix,sql,sqr] = scovint(x,0.95);
    end
    unc = max(abs(sqr - x0),abs(sql - x0));    
end