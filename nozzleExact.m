function [U_exact, RHO_exact, P_exact, M_exact] = nozzleExact(L, dx)    

%Constants
R=287;
g=1.4;
P01=100e3;
T01=300;
RHO01=P01/(R*T01);
C01=sqrt(g*P01/RHO01);
tol=1e-10;
cnt=0;
ss=0;%supersonic case (aka q2)
a=0;
xshock=7;

%Mesh info
xvar=dx;
xlength=L/xvar;
data=zeros(xlength,4);%column1-x, column2-S(x), column3=Mach, column4=Press

for x=0:xvar:L
    res=1;
    cnt=cnt+1;
    data(cnt,1)=x;
    Scrit=0.8;%q1 -- 0.8, q2 -- 1.0
    
    if (x<=(L/2) || ss==0)
        Mlo=0;
        Mhi=1;
        Mi=Mlo;
    elseif(x>(L/2) && ss==1)
        if(x<=xshock)
            Mlo=1;
            if (x==6.2)
                Mhi=1.25;
            else
                Mhi=200;
            end
            Mi=Mhi;
        elseif(x>xshock)
            Mlo=0;
            Mhi=1;
            Mi=Mlo;
        end
    end
    
    if (x<=(L/2))
        data(cnt,2)=1+1.5*(1-(data(cnt,1)/5))^2;
    elseif (x>(L/2))
        data(cnt,2)=1+0.5*(1-(data(cnt,1)/5))^2;
    end
    
    if (x>xshock && ss==1 && a==0)
        a=1;
        Prel_num=(((g+1)*data(cnt-1,3)^2/2)/(1+((g-1)*data(cnt-1,3)^2/2)))^(g/(g-1));
        Prel_den=((2*g*data(cnt-1,3)^2/(g+1))-((g-1)/(g+1)))^(1/(g-1));
        P0R=P01*Prel_num/Prel_den;
        rho01=P01/(R*T01);
        rho0R=P0R/(R*T01);
        a01=sqrt(g*P01/rho01);
        a0R=sqrt(g*P0R/rho0R);
        rhoL_aL_crit=rho01*a01*(2/(g+1))^((g+1)/(2*(g-1)));
        rhoR_aR_crit=rho0R*a0R*(2/(g+1))^((g+1)/(2*(g-1)));
        ScritR=Scrit*(rhoL_aL_crit/rhoR_aR_crit);
    end
    
    while(abs(res)>tol)
        if(x>xshock && ss==1)
            ls=data(cnt,2)/ScritR;
        else
            ls=data(cnt,2)/Scrit;
        end
        
        rs=(1/Mi)*((2/(g+1))*(1+((g-1)*Mi^2)/2))^((g+1)/(2*(g-1)));
        res=ls-rs;
                       
        if(x<=(L/2) || ss==0)
            if (res<0)
                Mlo=Mi;
            elseif (res>=0)
                Mhi=Mi;
            end
            Mi=Mlo+(Mhi-Mlo)/2;
        elseif(x>(L/2) && ss==1)
            if(x<=xshock)
                if (res>0)
                    Mlo=Mi;
                elseif (res<=0)
                    Mhi=Mi;
                end
                Mi=Mlo+(Mhi-Mlo)/2;
            elseif(x>xshock)
                if (res<0)
                    Mlo=Mi;
                elseif (res>=0)
                    Mhi=Mi;
                end
                Mi=Mlo+(Mhi-Mlo)/2;
            end
        end
    end
    
    data(cnt,3)=Mi;
    
    %Pressure Calculation
    if (x>xshock && ss==1)
        data(cnt,4)=P0R*(1+((g-1)*Mi^2)/2)^(-g/(g-1));
    else
        data(cnt,4)=P01*(1+((g-1)*Mi^2)/2)^(-g/(g-1));
    end
    
    %Other flow parameters
    data(cnt,5) = RHO01/(1+((g-1)/2)*data(cnt,3)^2)^(1/(g-1));%rho
    data(cnt,7) = sqrt((C01^2)/(1+((g-1)/2)*data(cnt,3)^2));%sound speed
    data(cnt,6) = data(cnt,3)*data(cnt,7);%flow velocity
    data(cnt,8) = T01/(1+((g-1)/2)*data(cnt,3)^2);%temp
    data(cnt,9) = log((data(cnt,4)/P01)*(RHO01/data(cnt,5))^g)*(g-1)^-1;%entropy
    
end

%Assign parameters to single vector
U_exact=data(:,6);
RHO_exact=data(:,5);
P_exact=data(:,4);
M_exact=data(:,3);

end



