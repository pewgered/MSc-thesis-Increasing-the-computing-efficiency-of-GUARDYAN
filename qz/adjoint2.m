tic
R=10;
n_r=21;
n_phi=19;
adj_Chi1=zeros(n_r,1);
var_Chi1=zeros(n_r,1);
adj_Psi1=zeros(n_r,n_phi);
var_Psi1=zeros(n_r,n_phi);
z0=(0:(n_r-1))/(n_r-1)*R;
phi0=(0:(n_phi-1))/(n_phi-1)*pi;
[sin(phi0)',cos(phi0)']
g=1;
for i=1:n_r
    [adj_Chi1(i),var_Chi1(i)]=S_nA_LE_in(z0(i),g,1e7);
    for j=1:n_phi
        [i,j]
        [adj_Psi1(i,j),var_Psi1(i,j)]=S_nA_LE_out(z0(i),phi0(j),g,1e6);
    end
end


adj_Chi2=zeros(n_r,1);
var_Chi2=zeros(n_r,1);
adj_Psi2=zeros(n_r,n_phi);
var_Psi2=zeros(n_r,n_phi);
g=2;
for i=1:n_r
    %[adj_Chi2(i),var_Chi2(i)]=S_nA_LE_in(z0(i),g,1e6);
    for j=1:n_phi
        [i,j]
        [adj_Psi2(i,j),var_Psi2(i,j)]=S_nA_LE_out(z0(i),phi0(j),g,1e6);
    end
end
%save('adj','adj_Psi2','var_Psi2')

toc
%plot(z0,adj_Chi)