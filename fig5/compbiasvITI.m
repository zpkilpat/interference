%
% compbiasvITI.m
%
% computes the bias as a function of the difference between the 
%

h = 0.1;    % threshold
tau = 100;   % time constant of plasticity 
b = 0.01;      % rate of plasiticty change
p0 = 3;     % plasticity target value
s = 0.005;

N = 2000;    % number of grid points
Nsim = 500;  % number of sims per point
ind = [1:N];
dt = 0.1;  % time step
T = 1300;    % end time
dx = 2*pi/N;    % space step
x = linspace(-pi,pi-dx,N)';  % space grid
nt = round(T/dt)+1; % number of time points
timey = linspace(0,T,nt);   % time vector
cmu = zeros(nt,1);
I=0; x2=0;
ITIs = linspace(100,500,21);
bias = zeros(21,1);
bhi = bias; blo = bias; bs = zeros(Nsim,1);
tvec = linspace(0,500,21);

x1 = pi/2;
for k=1:21, ITI=ITIs(k);
for l=1:Nsim, U = zeros(N,1); P = ones(N,1);
for j = 1:nt-1,
    I=0; tmod = j*dt;
	if 10<tmod & tmod<60, I = 2*cos(x-x1); end
    if 360<tmod & tmod<410, I = -5; end
    if 360+ITI<tmod & tmod<410+ITI, I = 2*cos(x-x2); end
    
    fu = heaviside(U-h);
    f1c = dx*P'.*cos(x')*fu;    f1s = dx*P'.*sin(x')*fu;
    nos = randn*cos(x)+randn*sin(x);
  
    Un = U + dt*(-U+f1c*cos(x)+f1s*sin(x)+I)+s*sqrt(dt)*nos;
    Pn = P + dt*(1-P+b*fu.*(p0-P))/tau;
    if j*dt==710+ITI, 
        [junk,mi] = max(Un); bs(l) = x(mi);
        break
    end
    U = Un; P = Pn; 
end
end
    bias(k)=sum(bs)/Nsim; bhi(k)=bias(k)+std(bs); blo(k)=bias(k)-std(bs);
end
degbias = 180*bias'/pi; deghi=180*bhi'/pi; deglo=180*blo'/pi; degvec=ITIs;
% hold on, plot(x1vec,bhi,'r','linewidth',2);
% hold on, plot(x1vec,blo,'r','linewidth',2);
clear stdpatch
vertx=[degvec'; degvec(end:-1:1)']; vertx(:,2)=[deghi'; deglo(end:-1:1)'];
stdpatch.Vertices=vertx;
stdpatch.Faces=[1:42];
stdpatch.FaceColor=[250 128 114]/255;
stdpatch.EdgeColor=[250 128 114]/255;
figure(1), hold on, patch(stdpatch);
plot(degvec,degbias,'r-o','markersize',12,'linewidth',4);
axis([100 500 0 15])
set(gca,'xtick',[0:100:500]);
set(gca,'fontsize',30);

% figure(2), hold on, plot(degvec,deghi-degbias,'r-o','markersize',12,'linewidth',4);
% set(gca,'xtick',[-180:90:180]);
% set(gca,'fontsize',30);


