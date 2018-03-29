

x = br_list;
y = s_dfx;


%th_init = [0.1 0.1 0];
%FF = inline(' th(1)*exp(xp*th(2)) ','xp','th');
%[fcomp,p,kvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr(x, y.^0.1, th_init, FF,0.0001,100);
%yf = FF(x,p).^10;

k = 20;
[pn,yf] = polyfit(x,y.^(1/k),2);
yf = 1.5*yf.yf.^k;


loglog(x,y)
hold on;
plot(x,yf,'r')
hold off;




