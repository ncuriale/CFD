function [R_bar] = smooth(R,nodes)

beta = 0.6;

% Smoothing coefficients
BD = zeros(3*nodes,3*nodes);
%  BD(1:3,3*nodes-2:3*nodes) = -beta*eye(3);
BD(1:3,1:3) = (1+2*beta)*eye(3);
BD(1:3,4:6) = -beta*eye(3);
for j=2:nodes-1
    BD(3*j-2:3*j,3*j-5:3*j-3) = -beta*eye(3);
    BD(3*j-2:3*j,3*j-2:3*j) = (1+2*beta)*eye(3);
    BD(3*j-2:3*j,3*j+1:3*j+3) = -beta*eye(3);
end
BD(3*nodes-2:3*nodes,3*nodes-5:3*nodes-3) = -beta*eye(3);
BD(3*nodes-2:3*nodes,3*nodes-2:3*nodes) = (1+2*beta)*eye(3);
%  BD(3*nodes-2:3*nodes,1:3) = -beta*eye(3);

% Smooth residual
R_bar = BD\R;













