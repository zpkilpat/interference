%
% complocaldistend.m
%
% computes the distribution of end points
%

h = 0.1;    % threshold
tau = 100;   % time constant of plasticity 
b = 0.01;      % rate of plasiticty change
p0 = 3;     % plasticity target value
s = 0.01;

N = 2000;    % number of grid points
Nsim = 5000;  % number of sims per point
ind = [1:N];
dt = 0.1;  % time step
T = 450;    % end time
dx = 2*pi/N;    % space step
x = linspace(-pi,pi-dx,N)';  % space grid
nt = round(T/dt)+1; % number of time points
timey = linspace(0,T,nt);   % time vector
I=0; x2=0;
endp = zeros(Nsim,1);
U = zeros(N,1); P = ones(N,1);

for l=1:Nsim, 
    x1=2*pi*rand-pi;    % uniform
%     x1=randn/100;           % peaked
%     x1=randn/25+pi/2;      % skewed
for j = 1:nt-1,
    I=0; tmod = j*dt;
    if tmod<50, I = -5; end
	if 150<tmod & tmod<200, I = 2*cos(x-x1); end
    
    fu = heaviside(U-h);
    f1c = dx*P'.*cos(x')*fu;    f1s = dx*P'.*sin(x')*fu;
    nos = randn*cos(x)+randn*sin(x);
  
    Un = U + dt*(-U+f1c*cos(x)+f1s*sin(x)+I)+s*sqrt(dt)*nos;
    Pn = P + dt*(1-P+b*fu.*(p0-P))/tau;
    U = Un; P = Pn;
    
end
    [junk,mi] = max(Un); endp(l) = x(mi)-x1;
end
degp = 180*endp/pi;
for j=1:5000, if degp(j)<-180, degp(j)=degp(j)+360; end, end
for j=1:5000, if degp(j)>180, degp(j)=degp(j)-360; end, end
[nx,xc]=hist(degp,50); dxc=xc(2)-xc(1); nx=nx/dxc/Nsim;
% bar(xc,nx,1);
% set(gca,'fontsize',28);
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
% colors = hsv(numel(nx));
colors = zeros(50,3);
colors(:,1)=1-abs(xc)/max(abs(xc));
colors(:,2)=0;
colors(:,3)=0.3;
for j = 1:numel(nx)
    bar(xc(j),nx(j),dxc,'parent', aHand, 'facecolor', colors(j,:),'edgecolor','none');
end
set(gca,'fontsize',28);
axis([-15 15 0 0.12])
set(gca,'ytick',[]);

exp = dxc*xc*nx'
xs = xc.^2;
vr = dxc*xs*nx'-exp^2;
stdd =  sqrt(vr)


