load("chi.mat")
load("psi.mat")

SigmaT=[0.2,2.0];
R=10;
n_r=100;
n_phi=30;

phi0=(0:n_phi-1)/n_phi*180;
for i=1:n_phi
    lgnd(i)=string(phi0(i));
end



QZ=zeros(2,n_r,n_phi);

figure
hold on

for i=1:n_phi
    QZ(1,:,i)=adj_Chi1./adj_Psi1(:,i);%*sigmaT(1);
    plot(z0,QZ(1,:,i))
end
legend(lgnd)
title('QZ1')

figure
set(gca,'YScale','log')
hold on
for i=1:n_phi
    plot(z0,QZ(1,:,i))
end
legend(lgnd)
title('QZ1')

figure
hold on
for i=1:n_phi
    plot(z0,max(QZ(1,:,i),1))
end
legend(lgnd)
title('QZ1')

figure
set(gca,'YScale','log')
hold on
for i=1:n_phi
    plot(z0,max(QZ(1,:,i),1))
end
legend(lgnd)
title('QZ1')


figure
hold on
for i=1:n_phi
    QZ(2,:,i)=adj_Chi2./adj_Psi2(:,i);%*sigmaT(2);
    plot(z0,QZ(2,:,i))
end
legend(lgnd)
title('QZ2')

figure
set(gca,'YScale','log')
hold on
for i=1:n_phi
    plot(z0,QZ(2,:,i))
end
legend(lgnd)
title('QZ2')

figure
hold on
for i=1:n_phi
    plot(z0,max(QZ(2,:,i),1))
end
legend(lgnd)
title('QZ2')

figure
set(gca,'YScale','log')
hold on
for i=1:n_phi
    plot(z0,max(QZ(2,:,i),1))
end
legend(lgnd)
title('QZ2')

figure
set(gca,'YScale','log')
hold on
for i=1:n_phi
    plot(z0,QZ(1,:,i))
end
for i=1:n_phi
    plot(z0,QZ(2,:,i))
end
title('QZ')

figure
set(gca,'YScale','log')
hold on
for i=1:n_phi
    plot(z0,QZ(1,:,i)*SigmaT(1));
end
for i=1:n_phi
    plot(z0,QZ(2,:,i)*SigmaT(2));
end
title('SigmaZ')

id_W1=0.0534./adj_Chi1;
id_W2=0.0534./adj_Chi2;

figure
set(gca,'YScale','log')
hold on
plot(z0,id_W1);
plot(z0,id_W2);
title('Ideal Weight')

figure
hold on
plot(z0,id_W1);
plot(z0,id_W2);
title('Ideal Weight')

figure
plot(z0,adj_Psi1(:,1),z0,adj_Chi1,z0,adj_Chi1./adj_Psi1(:,1))

figure
semilogy(z0,adj_Psi1(:,1),z0,adj_Chi1,z0,adj_Chi1./adj_Psi1(:,1))

%save('QZ','QZ')