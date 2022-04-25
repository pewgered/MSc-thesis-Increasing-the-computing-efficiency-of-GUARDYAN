load("psi.mat")

sigmaT=[0.2,2.0];

R=10;
n_r=100;
n_phi=30;
z0=(0:n_r-1)/n_r*R+R/n_r/2;
phi0=(0:n_phi-1)/n_phi*180+180/n_phi/2;

figure
hold on
for i=1:n_phi
    plot(z0,adj_Psi1(:,i))
    lgnd(i)=string(phi0(i));
end
legend(lgnd)
title('Psi1')
for i=1:9
    xline(i)
end

% figure
% set(gca,'YScale','log')
% hold on
% for i=1:n_phi
%     plot(z0,adj_Psi1(:,i))
% end
% legend(lgnd)
% title('Psi2')

figure
hold on
for i=1:n_phi
    plot(z0,adj_Psi2(:,i))
end
legend(lgnd)
title('Psi2')
for i=1:9
    xline(i)
end

% figure
% set(gca,'YScale','log')
% hold on
% for i=1:n_phi
%     plot(z0,adj_Psi2(:,i))
% end
% legend(lgnd)
% title('Psi2')
% 
% figure
% hold on
% for i=1:n_phi
%     plot(z0,adj_Psi1(:,i))
% end
% for i=1:n_phi
%     plot(z0,adj_Psi2(:,i))
% end
% title('Psi')
% 
% figure
% set(gca,'YScale','log')
% hold on
% for i=1:n_phi
%     plot(z0,adj_Psi1(:,i))
% end
% for i=1:n_phi
%     plot(z0,adj_Psi2(:,i))
% end
% 
% title('Psi')
% figure
% set(gca,'YScale','log')
% hold on
% for i=1:n_phi
%     plot(z0,var_Psi1(:,i))
% end
% title('var Psi1')
% 
% figure
% set(gca,'YScale','log')
% hold on
% for i=1:n_phi
%     plot(z0,var_Psi2(:,i))
% end
% title('var Psi2')

