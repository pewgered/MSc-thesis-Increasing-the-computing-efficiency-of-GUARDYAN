clc
clear all

tic

load("chi.mat")
load("psi.mat")
load("QZ.mat")

R=10;
n_r=100;
n_phi=30;
d_r=R/n_r;
d_phi=pi/n_phi;

SigmaT=[0.2,2.0];
p2a=0.4;
p12=0.55;
det1=7;
det2=8;
os=1.001;

Res=0;
Var=0;

n_sirs=1;

n_neutrons=1e4;
Res_hist=zeros(1,n_neutrons);
%XS
QS=max(QZ,1);
%QS=QZ;
QS=ones(2,n_r,n_phi);
for i=1:n_neutrons
    alive=true;
    inside=true;
    first=true;
    E_g=1;
    W=1;
    r_neutron=[0,0,0];
    r_n=0;
    n_coll=1;
    dir=getIsotrDir;
    xyz_i=floor(dir);
    while alive && inside
        r=-log(1-rand); % optical distance
        while r>0 && alive && inside
            if first
                first=false;
                r_i=1;
                phi_i=1;
            else
                % finding adjoint
                r_n=norm(r_neutron);
                if r_n>R
                    inside=false;
                    bigsad1=1;
                    break
                end
                phi=acos(r_neutron*dir'/r_n);
                phi_i=floor(phi/d_phi)+1;
                r_i=floor(r_n/d_r)+1;
            end
            %[dir,E_g,r_i,phi_i]
            SigmaM=QS(E_g,r_i,phi_i)*SigmaT(E_g);
            qZ=QZ(E_g,r_i,phi_i);
            qZ=1;
            Q=qZ*(1+SigmaT(E_g)/(SigmaM*os)*(1-qZ));
            Q=1;
            s=qZ/Q;
            
            % lattice
            dists=([xyz_i+[1,1,1],xyz_i]*d_r-[r_neutron,r_neutron])./[dir,dir];
            dist_sorted=sort(dists,'descend');
            lat_i=find(dists==dist_sorted(3));
            dist=dists(lat_i);
            
            % real path length
            xs=r/(SigmaM*os*s);
            if xs>=dist
                xyz_i(mod(lat_i-1,3)+1)=xyz_i(mod(lat_i-1,3)+1)+1-2*floor((lat_i-1)/3);
                
                r_neutron=r_neutron+dist*dir;
                r=r-dist*SigmaM*os*s;
                W=W*exp(-SigmaM*os*dist*(1-s));
            else
                r_neutron=r_neutron+xs*dir;
                r=0;
                W=W*exp(-SigmaM*os*xs*(1-s));
            end
            %[d_r*xyz_i,r_neutron]
        end
        W=W/s;
        
        % escape
        r_n=norm(r_neutron);
        if r_n>=R
            inside=false;
            bigsad2=1;
        end
        
        %         if E_g==2
        %             Res_hist(i)=r_n;
        %             break
        %         end
        
%         if n_coll==5
%             phi=acos(r_neutron*dir'/r_n);
%             phi_i=floor(phi/d_phi)+1;
%             r_i=floor(r_n/d_r)+1;
%             if E_g==1
%                 Res_hist(i)=W*adj_Psi1(r_i,phi_i);
%             elseif E_g==2
%                 Res_hist(i)=W*adj_Psi1(r_i,phi_i);
%             end
%             break
%         end
%         n_coll=n_coll+1;
        
        
        % mindig kisebb e mint 1?
        if inside && rand<Q*SigmaT(E_g)/(SigmaM*os)
            if E_g==2
                % detection
                if det1<r_n && r_n<det2
                    W=W/Q*p2a;
                    Res=Res+W;
                    Res_hist(i)=W;
                    Var=Var+W^2;
                    alive=false;
                else
                    W=W*(1-p2a);
                    %[E_g,r_n,W,SigmaM]
                end
            end
            % SIRS
            [E_g,dir,WW]=SIRS(E_g,r_neutron,adj_Psi1,adj_Psi2,n_sirs,d_r,d_phi,p12);
            W=W*WW;
        end
        W=W*(os*SigmaM-SigmaT(E_g))/(os*SigmaM-SigmaT(E_g)*Q);
        % roulette
        if W<0.001
            if rand<0.1
                W=W*10;
            else
                alive=false;
            end
        end
    end
end

Var=sqrt(Var/Res^2-1/n_neutrons);
Res=Res/n_neutrons;
[Res,Var]
%compare result to this
[adj_Chi1(1),var_Chi1(1)*sqrt(4e7/n_neutrons)]

nnz(Res_hist)
Res_hist=nonzeros(Res_hist);
figure
%histogram(Res_hist,51,'BinLimits',[0,10.2])
histogram(Res_hist)
toc

%batch
%kiszökések
%ideális súly adott helyen? W*adj=Res
