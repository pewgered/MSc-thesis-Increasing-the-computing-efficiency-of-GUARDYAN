tic
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
n_neutrons=1e6;
Res_hist=zeros(n_neutrons,1);
for i=1:n_neutrons
    E_group=1;
    r_neutron=[0,0,0];
    W=1.0;
    inside=true;
    alive=true;
    while inside && alive
        x=-1.0/sigmaT(E_group)*log(1-rand);  % free path
        r_neutron=r_neutron+x*getIsotrDir;  
        r_n=norm(r_neutron);                
        
        if r_n>=R                % escape
            inside=false;
            E_group=1;
            Esc=Esc+1;
            break
        end
        
        if E_group==1
            if rand<p12         % down scat
                E_group=2;
            end
            
        elseif E_group==2
            if r_n>=det1 && r_n<=det2   % detection
                if rand<=p2a            % abs
                    %W=W*p2a;
                    Res=Res+W;
                    Var=Var+W^2;
                    alive=false;
                    Res_hist(i)=W;
                end
            else
                W=W*(1-p2a);
            end
        end
        
        if W<=0.01 && alive
            if rand<=0.1                 % roulette
                W=W*10;
            else
                alive=false;
                Rou=Rou+1;
            end
        end
    end
end

    
Var=sqrt(Var/Res^2-1/n_neutrons);
Res=Res/n_neutrons;
[Res,Var]
[Esc,Rou]

nnz(Res_hist)
Res_hist=nonzeros(Res_hist);
% figure
% histogram(Res_hist)
% histogram(Res_hist,51,'BinLimits',[0,10.2])
toc