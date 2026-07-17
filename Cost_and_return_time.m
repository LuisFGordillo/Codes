clearvars
clf
warning off

a = 1; K = 1;
trho = 0.247;      % the centered value for the OU process rho (rho star)
lambda = 0.01; 
phi = 0.04;     % OU parameters

psi = 0.1; 
eps = 0.001;

P = 0.95;
k =-log(1-P);
T = 40000; N = 2^17; dt = T/N;
X0 = 0.432; %0.5*(1+sqrt(1-4*trho)); 
rho0 = 0.2435+rand*(0.255-0.243);%trho;                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial value of the process rho


r = 1; %%%%%%%%%%%%%%%%% Cost factor

xs = 0.5*(1+sqrt(1-4*trho));                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stable equilibrium
xu = 0.5*(1-sqrt(1-4*trho));                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stable equilibrium
J = [2*xs-3*xs^2-trho -xs; 0 -lambda];      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Jacobian aevaluated at xs
S = [psi*xs 0; 0 phi^2];                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% matrix sigma*sigma'

 %w = 0.05;                                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sensitivity matrix
        %SM = w*eye(2); 
        %SM = [[ 3.5833   -0.1766];[-0.1766    0.0125]];
 
     % Get the dimensions of matrix J
     n = size(J, 1);
     
     % This equation solves vec(J * W + W * J') = vec(-S)
     I = eye(n);
     lhs = kron(I,J) + kron(J,I);  % Left-hand side matrix for the linear system
     rhs = -S(:);                   % Right-hand side vectorized
     
     % Solve the linear system
     wtemp = lhs \ rhs;
     
     % Reshape the solution vector back into the matrix form
     SM = reshape(wtemp, n, n);
     
     SM(1,1) = 0.8;
     SM(2,2) = SM(1,1);
    
    
    %CM = -J-0.5*S*inv(SM);                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Feedback regulator matrix
    B = [1;0];
    BB = [1 0];
    P = eye(2)-B*BB; %[[0 0];[0,1]];
    
    
    CM = -0.5*BB*(S+J*SM+SM*J')*(P+eye(2))*inv(SM);

Cost = zeros(1,20);
TimetoEq = zeros(1,20);
NN = 300;                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Number of paths for each value of epsilon

IntCost = zeros(20,NN);
Time = zeros(20,NN);

for m = 1:20

    for i=1:NN

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Increments of Brownian motion
        dW = sqrt(dt)*randn(2,N);
        R = 4; Dt = R*dt; L = N/R;
    
    
    
        rho = zeros(1,L);                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ornstein-Uhlenbeck process
        rhotemp = rho0;
        X = zeros(1,L);                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Population process
        Xtemp = X0;
    
       
        
        
            Winc = zeros(2,L);
            for j=1:L    
                Winc(1,j) = sum(dW(1,R*(j-1)+1:R*j));
                Winc(2,j) = sum(dW(2,R*(j-1)+1:R*j));
            end
        
            for j = 1:L    
                Xtemp = Xtemp + Dt*(a*Xtemp^2*(1-Xtemp/K)-rhotemp*Xtemp)+eps*sqrt(psi*Xtemp)*Winc(1,j) + Dt*(CM(1,1)*(Xtemp-xs) + CM(1,2)*(rhotemp-trho));
                rhotemp = rhotemp - lambda*(rhotemp - trho)*Dt + eps*phi*Winc(2,j);
        
                X(j) = Xtemp;
                rho(j) = rhotemp;
            end
        
        rho = [rho0,rho];
        X=[X0,X];
        

        condition = abs(X - xs) < 0.005;
        k = find(condition, 1, 'first');

        if isempty(k)
            Y = X;  % If no element satisfies the condition, Y is the entire X
        else
            Y = X(1:k);  % Y includes elements up to and including the first satisfying element
        end
            rho = rho(1:k);


        IntCost(m,i) = real(sum((r*(CM(1)*(Y-xs)+CM(2)*(rho-trho)).^2)*Dt));
        Time(m,i) = k;

    end

    

    Cost(m) = mean(IntCost(m,:));
    
    TimetoEq(m) = mean(Time(m,:));

    eps = eps+0.001;

end

    save('IntCost.dat','IntCost');
    save ('Time.dat','Time');


figure(1)
hold on
plot(0.001*(1:20),Cost);
%plot(0.001*(1:20),St(1,:)+St(2,:),'r--');
%plot(0.001*(1:20),St(1,:)-St(2,:),'r--');
xlabel('\epsilon','FontSize',16)
ylabel('Cost','FontSize',16)
xticks([0 0.005 0.01 0.015 0.02])
title('(b)','FontSize',16)
legend({'$C$','$C\pm$ standard deviation'},'Interpreter','latex','Location','northwest','FontSize',12)

figure(2)
hold on
plot(0.001*(1:20),TimetoEq);
%plot(0.001*(1:20),TimetoEq(1,:)+TimetoEq(2,:),'r--');
%plot(0.001*(1:20),TimetoEq(1,:)-TimetoEq(2,:),'r--');
xlabel('\epsilon','FontSize',16)
ylabel('Return time','FontSize',16)
xticks([0 0.005 0.01 0.015 0.02])
title('(c)','FontSize',16)
legend({'$\mathrm{E}(\tau)$','$\mathrm{E}(\tau)\pm$ standard deviation'},'Interpreter','latex','Location','best','FontSize',12)