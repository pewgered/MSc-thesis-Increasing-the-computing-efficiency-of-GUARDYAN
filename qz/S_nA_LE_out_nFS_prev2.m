function [Res,Var] = S_nA_LE_out_nFS_prev2(a,b,phi1,phi2,g,n_neutrons,Chi1,var_Chi1,Chi2,var_Chi2,n_r) %,n_neutrons)

sigmaT=[0.2,2.0];
p2a=0.4;
p12=0.55;
R=10;
det1=7;
det2=8;

Res=0;
Var=0;
Esc=0;
Rou=0;
%n_neutrons=1e5;
Wdist=zeros(n_neutrons,1);
for i=1:n_neutrons
    E_group=g;
    z0=(a^3+rand*(b^3-a^3))^(1/3);              % sampling source point
    r_neutron=[0,0,z0];
    phi0=phi1+(phi2-phi1)*rand;                 % sampling direction
    dir=[sin(phi0),0,cos(phi0)];
    while 1
        x=-1.0/sigmaT(E_group)*log(1-rand);     % free path
        r_neutron=r_neutron+x*dir;
        r_n=norm(r_neutron);
        
        if r_n>=R                % escape
            Esc=Esc+1;
            break
        end
        
        if E_group==1
            if rand<p12
                E_g=2;
            else
                r_i=ceil(r_n/(R/n_r));
                Res=Res+Chi1(r_i);       %previous adjoint
                Var=Var+Chi1(r_i)^2;     %variation absolute error squared?
                break
            end
            
            
        elseif E_group==2
            r_i=ceil(r_n/(R/n_r));
            Res=Res+(1-p2a)*Chi2(r_i);       %previous adjoint
            Var=Var+((1-p2a)*Chi2(r_i))^2;   %variation absolute error squared?
            break
        end
        
      
        dir=getIsotrDir;

    end
end

Var=sqrt(Var/Res^2-1/n_neutrons);
Res=Res/n_neutrons;
%[Res,Var];
%[Esc,Rou];
end