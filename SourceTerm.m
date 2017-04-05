function [G] = SourceTerm(P, S_p, j)
% Calculates the value of the source term at node 'j'

G(1,1) = 0;
G(2,1) = P*S_p(j);
G(3,1) = 0;

    








