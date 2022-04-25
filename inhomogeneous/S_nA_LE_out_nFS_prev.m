function [Res,Var] = S_nA_LE_out_nFS_prev(a,b,phi1,phi2,g,n_neutrons,Chi1,var_Chi1,Chi2,var_Chi2,n_r) %,n_neutrons)

SigmaT=[0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1;
        2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0];
SigmaM=2.0;

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
    W=1.0;
    inside=true;
    alive=true;
    phi0=phi1+(phi2-phi1)*rand;                 % sampling direction
    dir=[sin(phi0),0,cos(phi0)];
    while inside && alive
        x=-1.0/SigmaM*log(1-rand);     % free path
        r_neutron=r_neutron+x*dir;
        r_n=norm(r_neutron);
        
        if r_n>=R                % escape
            inside=false;
            Esc=Esc+1;
            break
        end
        
        is=ceil(r_n);
        if rand<SigmaT(E_group,is)/SigmaM
            
            if E_group==1
                if rand<p12         % down scat
                    E_group=2;
                else
                    r_i=ceil(r_n/(R/n_r));
                    Res=Res+W*Chi1(r_i);       %previous adjoint
                    Var=Var+(W*Chi1(r_i))^2;   %variation absolute error squared?
                    alive=false;
                    break
                end
                
                
            elseif E_group==2
                if r_n>=det1 && r_n<=det2   % detection
                    if rand<=p2a            % abs
                        Res=Res+W;
                        Var=Var+W^2;
                        alive=false;
                        Wdist(i)=W;
                    end
                else
                    W=W*(1-p2a);
                    r_i=ceil(r_n/(R/n_r));
                    Res=Res+W*Chi2(r_i);       %previous adjoint
                    Var=Var+(W*Chi2(r_i))^2;   %variation absolute error squared?
                    alive=false;
                end
            end
        end
        
        if alive
            if W<=0.001
                if rand<=0.1                  % roulette
                    W=W*10;
                else
                    alive=false;
                    Rou=Rou+1;
                end
            end
        end
        
        if alive
            dir=getIsotrDir;
        end
    end
end

Var=sqrt(Var/Res^2-1/n_neutrons);
Res=Res/n_neutrons;
%[Res,Var];
%[Esc,Rou];
end