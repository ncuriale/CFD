function [res] = residual(Q, S, S_p, nodes, dx, gam, bcFlag)

for j=1:nodes
    
    [P, rho, u, M, T, c, e] = flowParam_node(S, Q, j);
    
    [G] = SourceTerm(P, S_p, j);
    [FD] = FluxDiff(j, Q, S, nodes, gam, dx, bcFlag);
    [D] = Dissipation(j, Q, S, gam, nodes, dx, bcFlag);
        
    % Calculate residual at node 'j' 
    res(3*j-2:3*j,1) = -(G - FD + D);
    
end

end

