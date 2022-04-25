load("chi.mat")
tic
R=10;
n_r=100;
% adj_Chi1=zeros(n_r,1);
% var_Chi1=zeros(n_r,1);
z0=(0:n_r-1)/n_r*R;
n_neutrons=1e8;
g=1;
% for i=1:n_r
%     [adj_Chi1(i),var_Chi1(i)]=S_nA_LE_in_nFS(z0(i),z0(i)+R/n_r,g,n_neutrons);
%     i
% end

%n_neutrons=1e7;
% adj_Chi2=zeros(n_r,1);
% var_Chi2=zeros(n_r,1);
g=2;
for i=21:35%100%60:90
    [adj_Chi2(i),var_Chi2(i)]=S_nA_LE_in_nFS(z0(i),z0(i)+R/n_r,g,n_neutrons);
    i
end

% for i=91:100
%     [adj_Chi2(i),var_Chi2(i)]=S_nA_LE_in_nFS_prev(z0(i),z0(i)+R/n_r,g,n_neutrons,adj_Chi2,var_Chi2,n_r);
%     %[adj_Chi2(i),var_Chi2(i)]
%     i
% end
% 
% 
% for i=59:-1:1
%     [adj_Chi2(i),var_Chi2(i)]=S_nA_LE_in_nFS_prev(z0(i),z0(i)+R/n_r,g,n_neutrons,adj_Chi2,var_Chi2,n_r);
%     %[adj_Chi2(i),var_Chi2(i)]
%     i
% end

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