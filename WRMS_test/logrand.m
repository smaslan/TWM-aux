function [rnd] = logrand(A_min,A_max,sz)
    rnd = 10.^(log10(A_min) + (log10(A_max) - log10(A_min))*rand(sz));
end