%load("chi.mat")
tic
R=10;
n_r=100;
adj_Chi1=zeros(n_r,1);
var_Chi1=zeros(n_r,1);
z0=(0:n_r-1)/n_r*R;

g=1;
n_neutrons=1e5;
parfor i=1:n_r
    [adj_Chi1(i),var_Chi1(i)]=S_nA_LE_in_nFS(z0(i),z0(i)+R/n_r,g,n_neutrons);
    i
end


adj_Chi2=zeros(n_r,1);
var_Chi2=zeros(n_r,1);
g=2;
%n_neutrons=1e4;
parfor i=1:n_r
    [adj_Chi2(i),var_Chi2(i)]=S_nA_LE_in_nFS(z0(i),z0(i)+R/n_r,g,n_neutrons);
    i
end
%save('chi','adj_Chi1','adj_Chi2','var_Chi1','var_Chi2')

toc
figure
plot(z0,adj_Chi1,z0,adj_Chi2)
figure
semilogy(z0,adj_Chi1,z0,adj_Chi2)
figure
plot(z0,var_Chi1,z0,var_Chi2)
figure
semilogy(z0,var_Chi1,z0,var_Chi2)