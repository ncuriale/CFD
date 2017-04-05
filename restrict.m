function [Q_2h, R_2h, S, S_p, nodes, dx] = restrict(Q_h, R_h, nodes, bcFlag)

nodes = (nodes-1)/2;
L = 10; 
dx = L/(nodes+1); 

% Calculate values of nozzle cross sectional area 
[S, S_p] = nozzleArea(nodes, dx, bcFlag);

% Q_h to Q_2h
Q_2h = zeros(3*length(nodes),1);
for j=1:nodes    
    Q_2h(3*j-2:3*j,1) = Q_h(6*j-2:6*j,1); 
end    

% R_h to R_2h
R_2h = zeros(3*length(nodes),1);
for j=1:nodes
    R_2h(3*j-2:3*j,1) = (1/4)*R_h(6*j-5:6*j-3,1) + (1/2)*R_h(6*j-2:6*j,1) + (1/4)*R_h(6*j+1:6*j+3,1);  
end 
    
   
    

