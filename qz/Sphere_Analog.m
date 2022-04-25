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

n_neutrons=1e6;
Res_hist=zeros(1,n_neutrons);
for i=1:n_neutrons
    E_group=1;
    r_neutron=[0,0,0];
    
    inside=true;
    alive=true;
    while inside && alive
        x=-1.0/sigmaT(E_group)*log(1-rand);  % free path
        r_neutron=r_neutron+x*getIsotrDir;
        r_n=norm(r_neutron);
        
%         if E_group==2
%             if r_n<R
%                 Res_hist(i)=r_n;
%             else
%                 Res_hist(i)=R+0.19;
%             end
%             break
%         end
        
        if r_n>=R                % escape
            inside=false;
            E_group=1;
            Esc=Esc+1;
            break
        end
        
        if E_group==1
            if rand<p12    % down scat
                E_group=2;
            end
            
        elseif E_group==2
            if rand<=p2a       % abs
                if r_n>=det1 && r_n<=det2
                    %Res_hist(i)=r_n;
                    Res=Res+1;
                    Var=Var+1^2;
                end
                alive=false;
            end
        end
    end
end
Var=sqrt(Var/Res^2-1/n_neutrons);
Res=Res/n_neutrons;
[Res,Var]
Esc

nnz(Res_hist)
Res_hist=nonzeros(Res_hist);
% figure
% histogram(Res_hist,51,'BinLimits',[0,10.2])
%histogram(Res_hist)


toc