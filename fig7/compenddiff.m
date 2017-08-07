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
T = 800;    % end time
dx = 2*pi/N;    % space step
x = linspace(-pi,pi-dx,N)';  % space grid
nt = round(T/dt)+1; % number of time points
timey = linspace(0,T,nt);   % time vector
I=0;
endp = zeros(Nsim,1);
U = zeros(N,1); P = ones(N,1);
nc = round((T-200)/dt)+1; st = round(200/dt);
cms = zeros(Nsim,nc);

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
    if 200<tmod,
        [junk,mi] = max(Un); cms(l,j+1-st) = x(mi)-x1;
        if cms(l,j+1-st)<-pi, cms(l,j+1-st)=cms(l,j+1-st)+2*pi; end
        if cms(l,j+1-st)>pi, cms(l,j+1-st)=cms(l,j+1-st)-2*pi; end
    end
end
end
dcen = 180*cms/pi;

cmv = var(dcen);
plot([0:dt:T-200],cmv,'k','linewidth',5);
set(gca,'xtick',[0:100:600]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);

