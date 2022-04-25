%clc
clear all
tic

load("chi.mat")
load("psi.mat")
load("QZ.mat")

%mean(adj_Psi1(1,:))

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
os=1.1;
%os=1;
%os=2;
%os=10;

Res=0;
Var=0;
Rou=0;
n_SIRS=10;
%n_SIRS=10;

resW=0;
varW=0;
nW=0;

n_neutrons=1e3;
Res_hist=zeros(1,n_neutrons);
W_1=zeros(1,n_neutrons);
W_1_z=zeros(1,n_neutrons);
W_2=zeros(1,n_neutrons);
W_2_z=zeros(1,n_neutrons);

W_3=zeros(1,n_neutrons);
W_3_z=zeros(1,n_neutrons);
W_4=zeros(1,n_neutrons);
W_4_z=zeros(1,n_neutrons);

%XS
QS=max(QZ,1);
%QS=QZ;
%QS=ones(2,n_r,n_phi);
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
    %dir=[sqrt(0.0000002) sqrt(0.0000001) sqrt(0.9999997)];
    xyz_i=floor(dir);
    while alive && inside
        r=-log(1-rand); % optical distance
        %r=10;
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
            SigmaM=QS(E_g,r_i,phi_i)*SigmaT(E_g);%*os; %os?
            qZ=QZ(E_g,r_i,phi_i);
            %qZ=1;
            %Q=qZ*(1+SigmaT(E_g)/(SigmaM*os)*(1-qZ));
            Q=qZ*SigmaM*os/(SigmaM*os+(qZ-1)*SigmaT(E_g));
            %Q=1.1;
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
        %             Res_hist(i)=n_coll;
        %             break
        %         end
        
        
        
        
        
        % mindig kisebb e mint 1?
        if inside && rand<Q*SigmaT(E_g)/(SigmaM*os)
            W=W/Q;%/os;
            
%             if n_coll==1
%                 r_n=norm(r_neutron);
%                 phi=acos(r_neutron*dir'/r_n);
%                 phi_i=floor(phi/d_phi)+1;
%                 r_i=floor(r_n/d_r)+1;
%                 if E_g==1
%                     %Res_hist(i)=r_n;
%                     W_1(i)=W;
%                     W_1_z(i)=r_n;
%                 elseif E_g==2
%                     %Res_hist(i)=r_n;
%                     W_2(i)=W;
%                     W_2_z(i)=r_n;
%                 end
%                 %break
%             end
            
%             if r_n<7
%                 resW=resW+W*adj_Chi1(r_i);
%                 varW=varW+(W*adj_Chi1(r_i))^2;
%                 nW=nW+1;
%             end
%             break
            

            if E_g==2
                % detection
                if r_n>det1 && r_n<det2
                    if rand<p2a
                        Res=Res+W;
                        Var=Var+W^2;
                        alive=false;
                        Res_hist(i)=W;
                        break
                    end
                else
                    W=W*(1-p2a);
                    %[E_g,r_n,W,SigmaM]
                end
            end
            
            % SIRS
            if n_SIRS~=1
                [E_g,dir,WW]=SIRS(E_g,r_neutron,adj_Psi1,adj_Psi2,n_SIRS,d_r,d_phi,p12);
                W=W*WW;
                %Res_hist(i)=WW;
                %if WW>1000
                %    [WW,r_n,n_coll]
                %end
            else
                if E_g==1
                    if rand<p12
                        E_g=2;
                    end
                end
                dir=getIsotrDir;
            end
            
            
%             if n_coll==1
%                 r_n=norm(r_neutron);
%                 phi=acos(r_neutron*dir'/r_n);
%                 phi_i=floor(phi/d_phi)+1;
%                 r_i=floor(r_n/d_r)+1;
%                 if E_g==1
%                     %Res_hist(i)=r_n;
%                     W_3(i)=W*adj_Psi1(r_i,phi_i);
%                     W_3_z(i)=r_n;
%                 elseif E_g==2
%                     %Res_hist(i)=r_n;
%                     W_4(i)=W*adj_Psi2(r_i,phi_i);
%                     W_4_z(i)=r_n;
%                 end
%                 break
%             end
            
            
            n_coll=n_coll+1;
        else
            %sulyfaktor=1/s
            W=W*(os*SigmaM-SigmaT(E_g))/(os*SigmaM-SigmaT(E_g)*Q);
        end
        
        
        %         if E_g==2
        %             Res_hist(i)=W;
        %             Res_hist(i)=r_n;
        %             %Res_hist(i)=n_coll;
        %             break
        %         end
        
        %         if W>1000
        %             [W,r_n]
        %         end
        
        
        
        %         roulette
%         if W<0.001
%             if rand<0.1
%                 W=W*10;
%             else
%                 alive=false;
%                 Rou=Rou+1;
%             end
%         end
        
        
%         n_coll=n_coll+1;
    end
end

% varW=sqrt(varW/resW^2-1/nW);
% resW=resW/nW;
% [resW,varW]
% 
% 
Var=sqrt(Var/Res^2-1/n_neutrons);
Res=Res/n_neutrons;
[Res,Var,Res*Var]
[0.0534,0.0012*sqrt(1e7/n_neutrons)]
Rou

% nnz(Res_hist)
% Res_hist=nonzeros(Res_hist);

% figure
% histogram(Res_hist)



%histogram(Res_hist,200,'BinLimits',[0,10])
% xline(7)
% xline(8)

% figure
% plot(W_1_z,W_1,'b.',W_2_z,W_2,'r.')

% id_W1=0.0534./adj_Chi1;
% id_W2=0.0534./adj_Chi2;
% z0=(0:n_r-1)/n_r*R+R/n_r/2;
% 
% 
% idealW=zeros(1,n_r);
% for i=1:n_r
%     sumqz=0;
%     for j=1:i
%         sumqz=sumqz+(1-QZ(1,j,1));  
%     end
%     idealW(i)=1/QZ(1,i,1)*exp(-SigmaT(1)*sumqz*d_r);
% end

% figure
% hold on
% set(gca,'YScale','log')
% plot(W_1_z,W_1,'b.',W_2_z,W_2,'r.')
% plot(z0,id_W1);
% plot(z0,id_W2);
% plot(z0,idealW);
% 
% figure
% hold on
% plot(W_1_z,W_1,'b.',W_2_z,W_2,'r.')
% plot(z0,id_W1);
% %plot(z0,id_W2);
% plot(z0,idealW);

% figure
% hold on
% plot(W_3_z,W_3,'b.',W_4_z,W_4,'r.')
% yline(0.0534)
% 
% figure
% hold on
% set(gca,'YScale','log')
% plot(W_3_z,W_3,'b.',W_4_z,W_4,'r.')
% yline(0.0534)

toc
%batch
%kiszökések
%ideális súly adott helyen? W*adj=Res
%pareto statisztika ellenőrzés
%chi reakció
%torzítatlanság?
%1,2,... ütközés után az adjungáltak átlaga (torzítva qz)
%1,2,... kcsh után kimenő adjungáltak átlaga
%1D 0 szórás

%eredmények diszkutálása
%qz csak s vagy csak Q vs os
%variance qZ, SIRS 
%mit lehet legjobbat kifozni belőle?
%homogén qz felbontás rontása
%nincs szögfüggés
%grid?
%futási idők
%inhomogén nagyon advanced