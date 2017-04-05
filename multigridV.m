function [Q_h] = multigridV(Q_h, S, S_p, nodes, dx, CFL, gam, bcFlag, tsFlag, num_mg, dt)

%Constants for iteration
q=5;% # of stages
N=1;% iterations of multistage 
    
% V cycle
% Return if coarsest grid
if (num_mg == 1)
    return
end

% Find residual based on restricted Q_h
R_h = residual(Q_h, S, S_p, nodes, dx, gam, bcFlag);
% Restrict Q and R to coarse grid
[Q_2h0, R_2h0, S, S_p, nodes, dx] = restrict(Q_h, R_h, nodes, bcFlag);
% Find residual based on restricted Q_2h0
res = residual(Q_2h0, S, S_p, nodes, dx, gam, bcFlag);
% Forcing function 'P_2h'
P_2h = R_2h0 - res;
% Complete # time steps to find updated Q_2h
Qi=Q_2h0;
for i=1:N
    Qi = iteration(Qi, S, S_p, nodes, dx, CFL, gam, bcFlag, tsFlag, P_2h, dt, q);
    if (i==N)
        Q_2h=Qi;
    end
end
% Perform next multigrid cycle
[Q_2h] = multigridV(Q_2h, S, S_p, nodes, dx, CFL, gam, bcFlag, tsFlag, num_mg-1, dt);
% Prolong error and update Q
Q_h = Q_h + prolong( (Q_2h - Q_2h0) , nodes);

end


                
                     