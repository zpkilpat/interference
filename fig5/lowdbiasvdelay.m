%
% lowdbiasvdelay.m
%
% simulates a bump in a single plastic layer with short term plasticity.
%

h = 0.1;    % threshold
tau = 100;   % time constant of plasticity 
b = 0.01;      % rate of plasiticty change
p0 = 3;     % plasticity target value
% a = pi/2-asin(h*(1+b)/(1+b*p0))/2;
a = pi/2-asin(h)/2;
C = b*(p0-1)/2/(1+b)/tan(a);

dt = 0.1;
T = 500;
nt = round(T/dt)+1;
tlow = linspace(0,T,nt);

oldamp = exp(-200/tau)*(1-exp(-(1+b)*350/tau));
oldamp = exp(-200/tau);
newamp = (1-exp(-(1+b)*50/tau));
xnew = 0;
xold = pi/2;
% Kold = sign(x-xold).*(1-cos(x-xold))-tan(a)*sin(x-xold);
% Knew = sign(x-xnew).*(1-cos(x-xnew))-tan(a)*sin(x-xnew);

y = zeros(nt,1); y(1) = xnew;
yp = zeros(nt,1); yp(1) = xnew;
yo = zeros(nt,1); yo(1) = xold;

for j=1:nt-1,
   y(j+1) = y(j)+dt*C*oldamp*heaviside(2*a-abs(y(j)-xold))*(sign(y(j)-xold).*(1-cos(y(j)-xold))-tan(a)*sin(y(j)-xold))*exp(-j*dt/tau)+...
       dt*C*(1-exp(-(1+b)*(50+j*dt)/tau))*heaviside(2*a-abs(y(j)-yp(j)))*(sign(y(j)-yp(j)).*(1-cos(y(j)-yp(j)))-tan(a)*sin(y(j)-yp(j)));
%    y(j+1) = y(j)+dt*C*oldamp*heaviside(2*a-abs(y(j)-xold))*(sign(y(j)-xold).*(1-cos(y(j)-xold))-tan(a)*sin(y(j)-xold))*exp(-j*dt/tau);
%    y(j+1) = y(j)+dt*C*oldamp*(sign(y(j)-xold).*(1-cos(y(j)-xold))-tan(a)*sin(y(j)-xold))*exp(-j*dt/tau)+...
%        dt*C*newamp*(sign(y(j)-yp(j)).*(1-cos(y(j)-yp(j)))-tan(a)*sin(y(j)-yp(j)));
   yp(j+1) = yp(j)-dt*(1+b)*(yp(j)-y(j))/tau;
%    yo(j+1) = yo(j)-dt*(1+b)*(yo(j)-y(j))/tau;
%    yo(j+1) = yo(j);
%    y(j+1) = y(j)+dt*C*oldamp*(sign(y(j)-xold).*(1-cos(y(j)-xold))-tan(a)*sin(y(j)-xold))*exp(-j*dt/tau);
end
lowdbias=y; dlowbias=180*lowdbias/pi;

figure(1), hold on, plot(tlow,dlowbias,'b','linewidth',4);