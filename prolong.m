function [E_h,nodes] = prolong(E_2h, nodes)

nodes = nodes*2 + 1;
L = 10; 
dx = L/(nodes+1); 

E_h = zeros(3*nodes,1);
for j=1:((nodes-1)/2)    
    E_h(6*j-2:6*j,1) = E_2h(3*j-2:3*j,1);    
end    

for j=1:((nodes-1)/2)+1  
    if j == 1 
        E_h(6*j-5:6*j-3,1) = 0.5*(E_2h(3*j-2:3*j,1));   
    elseif j <= ((nodes-1)/2)
        E_h(6*j-5:6*j-3,1) = 0.5*(E_2h(3*j-5:3*j-3,1) + E_2h(3*j-2:3*j,1));        
    elseif j == ((nodes-1)/2)+1    
        E_h(6*j-5:6*j-3,1) = 0.5*(E_2h(3*j-5:3*j-3,1));         
    end    
end  

    
   
    

