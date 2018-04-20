
function win = hanningw(M) % M: whole window, including start-end zeros
w = (1 - cos(2*pi*(0:(M-1))'/(M)));
%0:(M-1)
%Normalizing
win=(w*M/sum(w))';

