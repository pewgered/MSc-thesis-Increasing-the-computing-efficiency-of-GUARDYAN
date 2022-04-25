function [E_g,dir,W] = SIRS(E_g,r_neutron,adj_Psi1,adj_Psi2,n_sirs,d_r,d_phi,p12)

r_n=norm(r_neutron);
r_i=ceil(r_n/d_r);
if E_g==1
    % SIRS
    dirs=zeros(n_sirs,3);
    psis=zeros(n_sirs,1);
    E_gs=zeros(n_sirs,1);
    for ii=1:n_sirs
        dir=getIsotrDir;
        dirs(ii,:)=dir;
        phi=acos(r_neutron*dir'/r_n);
        phi_i=ceil(phi/d_phi);
        if rand<p12
            E_gs(ii)=2;
            psis(ii)=adj_Psi2(r_i,phi_i);
        else
            E_gs(ii)=1;
            psis(ii)=adj_Psi1(r_i,phi_i);
        end
    end
    sumpsi=sum(psis)
    psi_i=0;
    randpsi=rand*sumpsi;
    % chose one based on psi
    while randpsi>0
        psi_i=psi_i+1;
        randpsi=randpsi-psis(psi_i);
    end
    psi_i
    E_g=E_gs(psi_i);
    dir=dirs(psi_i,:);
    W=sumpsi/psis(psi_i)/n_sirs;
    
else
    % SIRS
    dirs=zeros(n_sirs,3);
    psis=zeros(n_sirs,1);
    for ii=1:n_sirs
        dir=getIsotrDir;
        dirs(ii,:)=dir;
        phi=acos(r_neutron*dir'/r_n);
        phi_i=ceil(phi/d_phi);
        psis(ii)=adj_Psi2(r_i,phi_i);
    end
    sumpsi=sum(psis);
    psi_i=0;
    randpsi=rand*sumpsi;
    % chose one based on psi
    while randpsi>0
        psi_i=psi_i+1;
        randpsi=randpsi-psis(psi_i);
    end
    dir=dirs(psi_i,:);
    W=sumpsi/psis(psi_i)/n_sirs;
    
end

