n=1e7;
a=1;
tic
for i=1:n
    r=rand(1,3);
    a=norm(r);
end
toc

tic
for i=1:n
    r=rand(1,3);
    a=sqrt(r(1)^2+r(2)^2+r(3)^2);
end
toc
tic
for i=1:n
    r=rand(1,3);
    a=r(1)^2+r(2)^2+r(3)^2;
end
toc

tic
for i=1:n
    r=rand(1,3);
    a=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3));
end
toc

tic
for i=1:n
    r=rand(1,3);
    a=r(1)*r(1)+r(2)*r(2)+r(3)*r(3);
end
toc