function [U_exact, RHO_exact, P_exact, M_exact] = shocktubeExact(L, t, dx) 

g=1.4;
alpha=(g+1)/(g-1);

pL=10^5;
rhoL=1;
aL=sqrt(g*pL/rhoL);
ML=0;

pR=10^4;
rhoR=0.125;
aR=sqrt(g*pR/rhoR);
MR=0;

tol=1e-10;
res=1;

Plo=0;
Phi=20;
P=Phi;
while (abs(res)>tol)
    ls=sqrt(2/(g*(g-1)))*(P-1)/sqrt(1+alpha*P);
    rs=(2/(g-1))*(aL/aR)*(1-(pR*P/pL)^((g-1)/(2*g)));
    res=ls-rs;
    
    if (res>0)
        Phi=P;
    elseif (res<=0)
        Plo=P;
    end
    P=Plo+(Phi-Plo)/2;
end

p2=P*pR;
rho2=rhoR*(1+alpha*P)/(alpha+P);
a2=sqrt(g*p2/rho2);

p3=p2;
rho3=rhoL*(p3/pL)^(1/g);
a3=sqrt(g*p3/rho3);

V=(2/(g-1))*aL*(1-(p3/pL)^((g-1)/(2*g)));
M2=V/a2;
M3=V/a3;
C=(P-1)*aR^2/(g*V);

x0=L/2;
xvar=dx;
xlength=L/xvar;
data=zeros(xlength,5);%column1-x column2-rho column3-M

cnt=0;
for x=0:xvar:L
    cnt=cnt+1;
    data(cnt,1)=x;
    
    if(x<=x0-aL*t)%left state
        data(cnt,2)=rhoL;
        data(cnt,3)=ML;
        
    elseif(x>x0-aL*t && x<=(x0+(V*(g+1)/2-aL)*t))%state 5
        u5=(2/(g+1))*((x-x0)/t + aL);
        a5=u5-(x-x0)/t;
        p5=pL*(a5/aL)^(2*g/(g-1));
        rho5=g*p5/a5^2;
        M5=u5/a5;
        data(cnt,2)=rho5;
        data(cnt,3)=M5;
        
    elseif(x>(x0+(V*(g+1)/2-aL)*t) && x<=(x0+V*t))%state 3
        data(cnt,2)=rho3;
        data(cnt,3)=M3;
        
    elseif(x>(x0+V*t) && x<=(x0+C*t))%state 2
        data(cnt,2)=rho2;
        data(cnt,3)=M2;
        
    elseif(x>(x0+C*t))%right state
        data(cnt,2)=rhoR;
        data(cnt,3)=MR;
    end
        
end
 
%Assign parameters to single vector
RHO_exact=data(:,2);
M_exact=data(:,3); 
P_exact=data(:,4);
U_exact=data(:,5);  
    
end


