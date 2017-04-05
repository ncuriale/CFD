% AER 1318
% Nathan Curiale
% 1002506781

clear all
clc
format long g

% Define constants
gam=1.4;% Specific heat ratio
R=287; % Specific gas constant for air
L=10; % Length of nozzle in meters
nodes=415; % Number of interior mesh points
dx=L/(nodes+1); % Node spacing
tol=1e-14; % Residual convergence tolerance
CFL =5; % Courant number

bcFlag = input('Specify the flow problem to be solved (1-2): ');
tsFlag = input('Specify constant time step or local time stepping (1-2): ');
mgFlag = input('Specify no multigrid or V or W multigrid cycle type (0-1-2): ');
if (mgFlag~=0)
    num_mg = input('Specify number of multi-grid levels (max 4 for W type): ');
else
    num_mg=0;
end

% Include boundary conditions
[uL, PL, rhoL, eL, SL, QL, cL, uR, PR, rhoR, eR, SR, QR, cR]= BCs(bcFlag);

% Calculate values of nozzle cross sectional area 
[S, S_p] = nozzleArea(nodes, dx, bcFlag);

% Calculate initial value of Q
Q = zeros(3*nodes,1);
for j=1:nodes
    Q(3*j-2) = rhoL*S(j);
    Q(3*j-1) = rhoL*uL*S(j);
    Q(3*j) = eL*S(j);
end
%  [Q] = IC(bcFlag);

% Flow parameters at intial Q
[P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);

%Time step type
if (tsFlag == 1) % Const time step
    maxspeed=0;
    for j=1:nodes
        if ((abs(u(j))+c(j))>maxspeed)
            maxspeed=abs(u(j))+c(j);
        end
    end
    dt = CFL*dx/maxspeed;
elseif (tsFlag == 2) % Local time step
    dt = 0;
end

% 5 Stage Runge-Kutta to find solution
iter=0;
tic %loop timer
for n=1:10000
    Q0 = Q;%Q before iteration function
    
    P_2h = zeros(3*nodes,1);
    
    % use 5-stage Runge-Kutta time marching
    q=5;
    [Q] = iteration(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, tsFlag, P_2h, dt, q);
    
    if (num_mg > 0) % Perform multigrid cycle
        if (mgFlag==1 || num_mg<=2)
            [Q] = multigridV(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, tsFlag, num_mg, dt);
        elseif (mgFlag==2 && num_mg==3)
            [Q] = multigridW3(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, tsFlag, num_mg, dt);
        elseif (mgFlag==2 && num_mg==4)
            [Q] = multigridW4(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, tsFlag, num_mg, dt);
        elseif (mgFlag==2 && num_mg>4)
            fprintf('ERROR: Enter less multi-grid levels for W cycle \n');
            break;
        end
    end
    dQ = Q - Q0;%Update dQ
    
    %store and check tolerence
    normRes(n) = norm(dQ);
    if (normRes(n)/normRes(1)<tol)
        break;
    end
    
    if (mod(n,10)==1)
        %show convergence
        fprintf('Iteration: %.0f', iter);
        fprintf('    Norm(res): %.14f', norm(dQ));
        fprintf('    normRes(n)/normRes(1): %.14f\n', normRes(n)/normRes(1));
%         [P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);
%         clf
%         plot(M)
%         drawnow
    end
    iter = iter+1;
    
end
toc %end loop timer
    
fprintf('Total Iterations: %4f\n', iter);

% Final flow solution quantities
[P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);

% Calculate exact solutions
if (bcFlag==1)   
    [U_exact, RHO_exact, P_exact, M_exact] = nozzleExact(L, dx);    
elseif (bcFlag == 2)    
    [U_exact, RHO_exact, P_exact, M_exact] = transnozzleExact(L, dx);
end

% Plot numerical and exact solutions
plotSol(P, M, rho, P_exact, M_exact, RHO_exact, bcFlag, nodes)

% Plot convergence
plotConv(normRes, bcFlag, iter)

% Error between exact and numerical solutions
M_error = zeros(nodes,1);
P_error = zeros(nodes,1);
for j=1:nodes
    M_error(j) = M_error(j) + abs(M(j)-M_exact(j+1));
    P_error(j) = P_error(j) + abs(P(j)-P_exact(j+1));
end
M_errorNorm = norm(M_error);
M_errorSum = sum(M_error);
P_errorNorm = norm(P_error);
P_errorSum = sum(P_error);
fprintf('Mach Error Norm: %4f\n', M_errorNorm);
fprintf('Mach Error Sum: %4f\n', M_errorSum);
fprintf('Pressure Error Norm: %4f\n', P_errorNorm);
fprintf('Pressure Error Sum: %4f\n', P_errorSum);
