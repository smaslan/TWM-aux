function [y]=sincm(x)
if x~=0 
    y=sin(x)/x;
else
    y=1;
end