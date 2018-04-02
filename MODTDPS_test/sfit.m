
mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% add variable variate lib: 
addpath([mfld filesep() 'var']);

load('data_no_harm.mat')

for k = 1:numel(res)
    res{k}.dpx = (res{k}.s_dpx.^2 + abs(res{k}.m_dpx).^2/3).^0.5;
    res{k}.dfx = (res{k}.s_dfx.^2 + abs(res{k}.m_dfx).^2/3).^0.5;
    res{k}.dAx = (res{k}.s_dAx.^2 + abs(res{k}.m_dAx).^2/3).^0.5;
    res{k}.dox = (res{k}.s_dox.^2 + abs(res{k}.m_dox).^2/3).^0.5;
end

r1 = var_get_results_mat(res,p,vr,'fs_rat','f0_per','bits',1);
r2 = var_get_results_mat(res,p,vr,'fs_rat','f0_per','bits',2);




return

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




