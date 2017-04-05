function [S, S_p] = nozzleArea(nodes, dx, bcFlag)
%   Defines cross-sectional area along tube depending on case selected.

S = zeros(nodes,1); % Holds nozzle cross sectional values

for j = 1:(5/dx)
    S(j) = 1 + 1.5*(1-(j*dx)/5)^2;
end
for j = ((5/dx)+1):nodes
    S(j) = 1 + 0.5*(1-(j*dx)/5)^2;
end

S_p = zeros(nodes,1); % Holds nozzle cross sectional derivatives 

for j = 1:(5/dx)
    S_p(j) = -3/5 + 3*j*dx/25;
end
for j = ((5/dx)+1):nodes
    S_p(j) = -1/5 + j*dx/25;
end
         
end


 