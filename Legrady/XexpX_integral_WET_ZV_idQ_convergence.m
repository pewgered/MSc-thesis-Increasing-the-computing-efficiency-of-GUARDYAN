NSamp=1E6;
a=10;
Sigma=1;
ResWET=0;
VarWET=0;
n_r=20;
os=1.0;
ExpW(1:n_r)=0;
VarW(1:n_r)=0;
NSample_per_bin(1:n_r)=0;
SigmaZ(1:n_r)=0;
SigmaS(1:n_r)=0;
x_co(1:n_r)=0;
for i=1:n_r
    x_coo=a/n_r*(i-1)+a/n_r/2; %to calculate the adjoint at interval middles
    x_co(i)=x_coo;
    SigmaZ(i)=x_coo*exp(-x_coo)/(exp(-x_coo)*(x_coo+1)-exp(-a)*(a+1));
    SigmaS(i)=max(Sigma,SigmaZ(i));
end
%plot(x_co,SigmaZ)
%figure
%plot(x_co,SigmaS)
for i=1:NSamp
    VirtCollision=1;
    Inside=1;
    W=1.0;    
    x=0;
    i_cell=1;
    while(VirtCollision && Inside)
        rr=rand;
        r=-log(1-rr); %optical distance
        inTravel=1;
        while(inTravel && i_cell<=n_r)
            SigmaM=SigmaS(i_cell);
            % to set the ideal paramters
            qZ=SigmaZ(i_cell)/Sigma;
            Q=qZ*(1+Sigma/(SigmaM*os)*(1-qZ)); % Now this is the new stuff!
            s=qZ/Q; % this ensures that the ideal qz is sampled
            
            %s should also depend on i_cell later
            xs=r/(SigmaM*os*s); %
            %ez is a lényeg
            if(x+xs>=a/n_r*i_cell)
                step=a/n_r*i_cell-x;
                x=x+step;
                W=W*exp(-SigmaM*os*step*(1-s)); %
                r=r-step*SigmaM*os*s;%
                i_cell=i_cell+1;
            else
                x=x+xs;
                W=W*exp(-SigmaM*os*xs*(1-s));%
                inTravel=0;
            end
        end
        W=W/s;
        if(x>=a)
            Inside=0;
        end
        if(x<a && rand<Q*Sigma/(SigmaM*os))
            VirtCollision=0;
            W=W/Q;
            Contribution=x*W;
            ResWET=ResWET+Contribution;
            VarWET=VarWET+Contribution^2;
            % tally W distribution
            %ii=ceil(x/a*n_r);
            %ExpW(ii)=ExpW(ii)+W;
            %VarW(ii)=VarW(ii)+W^2;
            %NSample_per_bin(ii)=NSample_per_bin(ii)+1;
        end
        W=W*(os*SigmaM-Sigma)/(os*SigmaM-Sigma*Q);
        if W~=0
            jonapot=1
        end
    end
end
VarW=sqrt(VarW./ExpW.^2-1./NSample_per_bin);
ExpW=ExpW./NSample_per_bin;
VarWET=sqrt(VarWET/ResWET^2-1/NSamp);
ResWET=ResWET/NSamp;
ResIdeal=1-exp(-a)*(a+1);
fprintf('Ideal= %f WET= %f relvarWET= %f\n', ResIdeal,ResWET,VarWET);
% 
% figure;
% xx=1:n_r;
% plot(xx,ExpW,xx,exp(-Sigma*xx/n_r*a*(1-qZ))/qZ);
% title('W expectation - the ideal W');
% 
% figure;
% plot(xx,VarW);
% title('Relative Variance of Ws');