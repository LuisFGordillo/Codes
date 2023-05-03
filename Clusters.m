%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code generates the figures showing
% the emergence of vegetation clusters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L=15;           % lenght of the interval = 2*L
np=1500;        % label points in space
h=2*L/np;
x=linspace(-L,L-h,np);
T=6000;         % total time

D=0.02;         % death rate   
K=10;
dt=1/10;
I=0.055;
p=0.14;         % probability of high precipitation  


%**************************************
% facilitation-inhibition influence kernel
%**************************************

%a=2; b=3; c=1; d=1;
a1=1.1; b1=9; c1=1.05; d1=8;    % low
%a2=1.5; b2=9; c2=1.05; d2=7;    % high 1
a2=1.3; b2=9; c2=1.05; d2=8;       % high 2

PH1=(a1*exp(-b1*x.^2)-c1*exp(-d1*x.^2));
PH11=abs(PH1)>=0.0001;
PH=(PH1).*PH11;

PH2=(a2*exp(-b2*x.^2)-c2*exp(-d2*x.^2));
PH22=abs(PH2)>=0.0001;
PHr=(PH2).*PH22;

B=[PH(np/2+1:np),PH(1:np/2)];
Br=[PHr(np/2+1:np),PHr(1:np/2)];


%**************************************
% initial conditions
%**************************************

u0=5*(exp(-0.01*x.^2))';
%u0=10*ones(length(x),1);
%u0=4*rand(length(x),1);
%u0=zeros(length(x),1);
%u0=(1+sqrt(1-2*D/(K*I)))*K/2*ones(length(x),1)+0.1*rand(length(x),1);


s0=u0>=0.001;
S0=u0.*s0;
U=S0;
C=circulant(B);
Cr=circulant(Br);

for i=1:T
    
    if rand<=p
        Ctemp=Cr;
    else
        Ctemp=C;
    end
        
    R=Ctemp*S0*h;
    Stemp=S0+(R.*S0.*(1-S0/K)-D*S0)*dt;
    Rtemp=Ctemp*Stemp*h;
    Stemp2=(Rtemp.*Stemp.*(1-Stemp/K)-D*Stemp);
    Stemp1=(Stemp-S0)/dt;
    S1=S0+(dt/2)*(Stemp1+Stemp2);
    K1=S1>=0.0001;
    u=K1.*S1;
    U=[U u];
    S0=u;
end

[X,Y] = meshgrid(x,(1:T+1)*dt);
surf(X,Y,U')
axis tight
shading flat
view(0,90)
colorbar
caxis([0 10])
xlabel('space')
ylabel('time')