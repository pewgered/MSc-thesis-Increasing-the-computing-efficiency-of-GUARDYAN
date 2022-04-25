n=100;
a=rand(n,1);
vara=a.^2;
suma=sum(a);
vara=sqrt(sum(vara)/suma^2-1/n);
suma=suma/n;
[suma,vara]

m=400;
b=rand(m,1);
varb=b.^2;
sumb=sum(b);
varb=sqrt(sum(varb)/sumb^2-1/m);
sumb=sumb/m;
[sumb,varb]

k=1000;
sumc=0;
varc=0;
for i=1:k
    if rand<0.4
        sumc=sumc+suma;
        varc=varc+suma^2;
    else
        sumc=sumc+sumb;
        varc=varc+sumb^2;
    end
end
varc=sqrt(sum(varc)/sumc^2-1/k);
sumc=sumc/k;
[sumc,varc]