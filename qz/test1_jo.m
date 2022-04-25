tic
SigmaT=0.2;
SigmaM=0.5;
Q=0.6;
s=0.7;
os=3;
dx=0.1;
n=200;
res_hist=zeros(n,1);
res_n=zeros(n,1);
n_neutron=1e6;
for i=1:n_neutron
    r_n=0;
    W=1;
    virtual=true;
    while virtual
        r=-log(1-rand);
%                 xs=r/(SigmaM*s*os);
%                 while r>0
%                     if xs>=dx
%                         r_n=r_n+dx;
%                         r=r-dx*SigmaM*s*os;
%                         W=W*exp(-SigmaM*os*dx*(1-s));
%                     else
%                         r_n=r_n+xs;
%                         r=0;
%                         W=W*exp(-SigmaM*os*xs*(1-s));
%                     end
%                 end
        
%         W=W*exp(-r/s*(1-s));
        x=r/(SigmaM*os*s);
        W=W*exp(-x*SigmaM*os*(1-s));
        r_n=r_n+x;
        
        W=W/s;
        if rand<Q*SigmaT/(SigmaM*os)
            virtual=false;
            W=W/Q;
            ii=ceil(r_n/dx);
            if ii<n
                res_hist(ii)=res_hist(ii)+W;
                res_n(ii)=res_n(ii)+1;
            end
        else
            W=W*(os*SigmaM-SigmaT)/(os*SigmaM-SigmaT*Q);
        end
    end
end

%res_n'
%res_hist'
x=dx:dx:n*dx;
y=SigmaT*exp(-SigmaT*x);
figure
semilogy(x,res_hist/n_neutron/dx,x,y)

y=SigmaT*s*Q*exp(-SigmaT*s*Q*x);
figure
plot(x,res_n/n_neutron/dx,x,y)
toc