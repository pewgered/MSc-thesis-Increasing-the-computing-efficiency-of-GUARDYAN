function [Res,Var] = S_nA_LE_in_nFS_prev(a,b,g,n_neutrons,Chi,var_Chi,n_r) %,n_neutrons)

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
%n_neutrons=1e6;
Wdist=zeros(n_neutrons,1);
for i=1:n_neutrons
    E_group=g;
    z0=(a^3+rand*(b^3-a^3))^(1/3);
    r_neutron=[0,0,z0];
    r_n=z0;
    W=1.0;
    inside=true;
    alive=true;
    while inside && alive
        
        if E_group==1
            if rand<p12         % down scat
                E_group=2;
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
                if Chi(r_i)~=0
                    Res=Res+W*Chi(r_i);   %previous adjoint
                    Var=Var+(W*Chi(r_i))^2;   %variation absolute error squared?
                    alive=false;
                end
            end
        end
        
        if W<=0.1 && alive
            if rand<=W                  % roulette
                W=1.0;
            else
                alive=false;
                Rou=Rou+1;
            end
        end
        
        x=-1.0/sigmaT(E_group)*log(1-rand);  % free path
        r_neutron=r_neutron+x*getIsotrDir;  
        r_n=norm(r_neutron);                
        
        if r_n>=R                % escape
            inside=false;
            Esc=Esc+1;
            break
        end
        
    end
end

Var=sqrt(Var/Res^2-1/n_neutrons);
Res=Res/n_neutrons;
%[Res,Var];
%[Esc,Rou];
end