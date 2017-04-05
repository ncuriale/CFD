function [Q] = iteration(Q, S, S_p, nodes, dx, CFL, gam, bcFlag, tsFlag, P_2h, dt, q)

% Coefficients of each stage of the Runge-Kutta time marching scheme
alpha(1) = 1/4; 
alpha(2) = 1/6;
alpha(3) = 3/8; 
alpha(4) = 1/2;
alpha(5) = 1;
                                           
% Artificial dissipation weights
Dw = zeros(5,5);
Dw(1,1) = 1;
Dw(2,1) = 1;
Dw(3,1) = (1-0.56);
Dw(3,3) = 0.56;
Dw(4,1) = (1-0.56);
Dw(4,3) = 0.56;
Dw(5,1) = (1-0.56)*(1-0.44);
Dw(5,3) = 0.56*(1-0.44);
Dw(5,5) = 0.44;
                                           
Qi = Q; % Q at the current timestep    

for m=1:q % Time-marching stages
    for j=1:nodes
        % Flow parameters at current node
        [P, rho, u, M, T, c, e] = flowParam_node(S, Q, j);
        
        [G] = SourceTerm(P, S_p, j);% source term
        [FD] = FluxDiff(j, Q, S, nodes, gam, dx, bcFlag);% flux diff vector
        [D] = Dissipation(j, Q, S, gam, nodes, dx, bcFlag);% artificial dissipation
         
        % Artificial dissipation terms at stage 'm'
        Dsum = zeros(3,1);
        for n=1:m
            Dsum = Dsum + Dw(m,n)*D;
        end
        
        % Local timestep 
        if (tsFlag == 2) 
            dt = CFL*dx/(abs(u)+c);
        end
        
        res(3*j-2:3*j,1) = -(G - FD + Dsum);
        R(3*j-2:3*j,1) = alpha(m)*dt*(res(3*j-2:3*j,1) + P_2h(3*j-2:3*j,1));
  
    end
    
    % Residual smoothing
    R = smooth(R, nodes);
    
    % Update Q
    Q = Qi - R;

end

end

