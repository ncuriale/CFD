function[uL, PL, rhoL, eL, SL, QL, cL, uR, PR, rhoR, eR, SR, QR, cR] = BCs(bcFlag)

gam = 1.4; 

if (bcFlag == 1) 

    uL = 65.4510; 
    PL = 97534.31;
    rhoL = 1.1409;
    eL = PL/(gam-1) + 0.5*rhoL*uL^2;
    cL = sqrt(gam*PL/rhoL); 
    SL = 2.5; 
    QL =[rhoL*SL; rhoL*uL*SL; eL*SL];
    uR = 113.0563;
    PR = 92772.07; 
    rhoR = 1.1008; 
    eR = PR/(gam-1) + 0.5*rhoR*uR^2; 
    cR = sqrt(gam*PR/rhoR);
    SR = 1.5;
    QR =[rhoR*SR; rhoR*uR*SR; eR*SR]; 
    
elseif (bcFlag == 2)
    
    uL = 82.6934; 
    PL = 96084.92;
    rhoL = 1.1287;
    eL = PL/(gam-1) + 0.5*rhoL*uL^2;  
    cL = sqrt(gam*PL/rhoL);
    SL = 2.5; 
    QL =[rhoL*SL; rhoL*uL*SL; eL*SL]; 
    uR = 151.6196; 
    PR = 84973.94; 
    rhoR = 1.0538; 
    eR = PR/(gam-1) + 0.5*uR^2;  
    cR = sqrt(gam*PR/rhoR); 
    SR = 1.5; 
    QR =[rhoR*SR; rhoR*uR*SR; eR*SR];  

end


