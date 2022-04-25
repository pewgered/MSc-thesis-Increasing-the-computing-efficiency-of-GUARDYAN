tic
load("chi.mat")

R=10;
n_r=100;
n_phi=30;
z0=(0:n_r-1)/n_r*R;
phi0=(0:n_phi-1)/n_phi*pi;

% adj_Psi1=zeros(n_r,n_phi);
% var_Psi1=zeros(n_r,n_phi);
n_neutrons=1e6;
% g=1;
% parfor i=1:n_r
%     for j=1:n_phi
%         [adj_Psi1(i,j),var_Psi1(i,j)]=S_nA_LE_out_nFS_prev(z0(i),z0(i)+R/n_r,phi0(j),phi0(j)+pi/n_phi,g,n_neutrons,adj_Chi1,var_Chi1,adj_Chi2,var_Chi2,n_r);
%     end
%     i
% end

adj_Psi2=zeros(n_r,n_phi);
var_Psi2=zeros(n_r,n_phi);
%n_neutrons=1e4;
g=2;
parfor i=1:n_r
    for j=1:n_phi
        [adj_Psi2(i,j),var_Psi2(i,j)]=S_nA_LE_out_nFS_prev(z0(i),z0(i)+R/n_r,phi0(j),phi0(j)+pi/n_phi,g,n_neutrons,adj_Chi1,var_Chi1,adj_Chi2,var_Chi2,n_r);
    end
    i
end
%save('psi','adj_Psi1','adj_Psi2','var_Psi1','var_Psi2')

toc

figure
plot(z0,adj_Psi1(:,1))
title('Psi1')
