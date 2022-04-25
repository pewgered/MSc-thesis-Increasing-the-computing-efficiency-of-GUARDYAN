function [v] = getIsotrDir()

v=zeros(1,3);
rho2=1;

while rho2 >=1
    a=2*rand()-1;
    b=2*rand()-1;
    rho2=a^2+b^2;      
end
r=sqrt(1-rho2);

v(1)=1-2*rho2;
v(2)=2*a*r;
v(3)=2*b*r;

end