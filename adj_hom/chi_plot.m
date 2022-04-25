load("chi.mat")

R=10;
n_r=100;
z0=(0:n_r-1)/n_r*R+R/n_r/2;

figure
plot(z0,adj_Chi1,z0,adj_Chi2)
figure
semilogy(z0,adj_Chi1,z0,adj_Chi2)
figure
plot(z0,var_Chi1,z0,var_Chi2)
figure
semilogy(z0,var_Chi1,z0,var_Chi2)