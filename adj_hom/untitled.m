function [] = untitled()
n=1e7;
tic
for i=1:n
    x=sin(2*pi*rand);
    y=sqrt(1-x*x);
end
toc
tic
v=zeros(1,3);
for i=1:n
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
toc
tic
for i=1:n
    getIsotrDir;
end
toc
tic
for i=1:n
    getIsotrDir2D;
end
toc
end

