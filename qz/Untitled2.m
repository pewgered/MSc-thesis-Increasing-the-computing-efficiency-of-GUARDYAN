
vars=zeros(1,20);
parfor i=1:20
    vars(i)=(qz_sphere_test2(i,2))^2;
    i
end

vars1=zeros(1,20);
parfor i=1:20
    vars1(i)=(qz_sphere_test2(i,1))^2;
    i
end

figure
hold on 
plot(1:20,vars1,'LineWidth',2)
plot(1:20,vars,'LineWidth',2)
yline(0.0042^2,'LineWidth',2)
xlabel('Number of SIRS samples')
ylabel('Relative Variance')
ylim([0,2e-5])
legend('without qZ','with qZ','analog')