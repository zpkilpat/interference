%
% lowdvar_nofac.m
%
% simulates a bump in a single plastic layer with short term plasticity.
%

h = 0.1;    % threshold
s=0.01;
na = s/sqrt(2)/sqrt(1+sqrt(1-h^2));

N=5000;

dt = 0.1;
T = 600;
nt = round(T/dt)+1;
timey = linspace(0,T,nt);

cmm = zeros(N,nt);
xnew = 0;

for k=1:N, y = zeros(nt,1); 
for j=1:nt-1,
   y(j+1) = y(j)+sqrt(dt)*na*randn;
end
    cmm(k,:) = y;
end
degcm=180*cmm/pi;

cmv = var(degcm);
plot([0:dt:T],cmv,'b','linewidth',8);
set(gca,'xtick',[0:100:600]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);