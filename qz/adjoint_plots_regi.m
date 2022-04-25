load('adj.mat')

sigmaT=[0.2,2.0];

R=10;
n_r=21;
n_phi=19;
z0=(0:(n_r-1))/(n_r-1)*R;
phi0=(0:(n_phi-1))/(n_phi-1)*pi/2;

%for i=1:n_phi
%    adj_Psi1(1,i)=adj_Chi1(1);
%    adj_Psi2(1,i)=adj_Chi2(1);
%end

figure(1)
plot(z0,adj_Chi1)
title('chi1')

figure(2)
plot(z0,adj_Chi2)
title('chi2')

figure(3)
semilogy(z0,adj_Chi2)
title('chi2')

figure(4)
hold on
for i=1:n_phi
    plot(z0,adj_Psi1(:,i))
    lgnd(i)=string((i-1)/(n_phi-1)*180);
end
legend(lgnd)
title('Psi1')

figure(5)
hold on
for i=1:n_phi
    plot(z0,adj_Psi2(:,i))
    lgnd(i)=string((i-1)/(n_phi-1)*180);
end
legend(lgnd)
title('Psi2')

figure(6)
set(gca,'YScale','log')
hold on
for i=1:n_phi
    semilogy(z0,adj_Psi2(:,i))
    lgnd(i)=string((i-1)/(n_phi-1)*180);
end
legend(lgnd)
title('Psi2')

figure(7)
hold on
SigmaZ1=zeros(n_r,n_phi);
for i=1:n_phi
    SigmaZ1(:,i)=adj_Chi1./adj_Psi1(:,i);%*sigmaT(1);
    plot(z0,SigmaZ1(:,i))
end
legend(lgnd)
title('SigmaZ1')

figure(8)
hold on
SigmaZ2=zeros(n_r,n_phi);
for i=1:n_phi
    SigmaZ2(:,i)=adj_Chi2./adj_Psi2(:,i);%*sigmaT(2);
    plot(z0,SigmaZ2(:,i))
end
legend(lgnd)
title('SigmaZ2')

figure(9)
set(gca,'YScale','log')
hold on
SigmaZ2=zeros(n_r,n_phi);
for i=1:n_phi
    SigmaZ2(:,i)=adj_Chi2./adj_Psi2(:,i);%*sigmaT(2);
    plot(z0,SigmaZ2(:,i))
end
legend(lgnd)
title('SigmaZ2')

save('SigmaZ','SigmaZ1','SigmaZ2')