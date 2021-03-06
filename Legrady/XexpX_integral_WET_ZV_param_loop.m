clc
clear all
close all
% integral of versions of x exp(-x) as I lost the original
% Woodcock exponential transform test
NSamp=1E6;
a=10;
Diff=0.2;
ResWET=0;
VarWET=0;
NWdiag=100;
N_Adj_Res=20;
N_qs=200;
N_o=5;
SigmaX(1:N_Adj_Res+1)=1.0;
SigmaZ(1:N_Adj_Res+1)=0;
SigmaS(1:N_Adj_Res+1)=0;
for i=1:N_Adj_Res+1
    SigmaZ(i)=1.0; %to be filled
    %now let us make things nasty! Let us have varying majorant
    if (mod(i,2)==0)
        SigmaS(i)=2.2; % to be set
    else
        SigmaS(i)=1.2;
    end
end

Sigma=1.0;

i_HalfSigma=ceil(Sigma/a*0.5*NWdiag);
Sc_H=[];
i_OneSigma=ceil(Sigma/a*1.0*NWdiag);
Sc_1=[];
i_TwoSigma=ceil(Sigma/a*2.0*NWdiag);
Sc_2=[];
i_FiveSigma=ceil(Sigma/a*5.0*NWdiag);
Sc_5=[];
ExpWeight(1:NWdiag)=0;
VarWeight(1:NWdiag)=0;
%ResIdeal=p1*(1-exp(-a*XSigma(1))*(a*XSigma(1)+1))/XSigma(1)^2+(1-p1)*(1-exp(-a*XSigma(2))*(a*XSigma(2)+1))/XSigma(2)^2;
%ResIdeal=ResIdeal/a;
ResIdeal=1-exp(-a)*(a+1);
NSample_per_bin(1:NWdiag)=0;
% expected value and distributiontest for the Woodcock Exponential
% Transform
qZ=2.0; % the ideal multiplication foactor to Sigma_local Sigma_z=Sigma(x)*qz

q_a(N_o,1:N_qs)=0;
s_a(N_o,1:N_qs)=0;
Res_One(1:N_o,1:N_qs)=0;
Var_One(1:N_o,1:N_qs)=0;
Ns_One(1:N_o,1:N_qs)=0;
Res_Five(1:N_o,1:N_qs)=0;
Var_Five(1:N_o,1:N_qs)=0;
Ns_Five(1:N_o,1:N_qs)=0;

for i_o=1:N_o
    for i_qs=1:N_qs
        OverSample=i_o; % oversampling is set in a fashion that it does not affect the expected value of the cross section sampled at the end
        % virtual/real collision selsection chance
        % max possible q to stay below 1 is:
        SigmaMaj=SigmaS(1); % this is supposed to account for assuming a lack of a tight majorant at every phase space loation
% special spots for histogram. Any further modification to this value
% happens through the oversampling parameter. SigmaMaj=Sigma would make Exp
% Transformation alone the ideal choice.
        q_max=SigmaMaj*OverSample/Sigma;
        q_max=1.8; % now zooming in
        q_min=0.05; % just not to make it 0
        q_min=1.0; % now zooming in
        q=(q_max-q_min)/N_qs*i_qs+q_min;
        if(q*Sigma/(SigmaMaj*OverSample)>1)
            q=1.0;
        end
        q_a(i_o,i_qs)=q;
        %q=qZ;
        s=qZ/q; % this ensures that the ideal qz is sampled
        s_a(i_o,i_qs)=s;
        % exponential transform parameter
        fprintf('qz = %f o= %f q= %f s= %f\n', qZ,OverSample,q,s);
        
        for i=1:NSamp
            VirtCollision=1;
            weight=1.0;
            Inside=1;
            x=0;
            i_cell=1;
            while(VirtCollision && Inside)
                r=rand;
                r=-log(1-r); %optical distance
                inTravel=1;
                while(inTravel && i_cell<=N_Adj_Res)
                    SigmaMaj=SigmaS(i_cell);
                    %s should also depend on i_cell later
                    xs=r/(SigmaMaj*s*OverSample);
                    
                    if(x+xs>=a/N_Adj_Res*i_cell)
                        step=a/N_Adj_Res*i_cell-x;
                        x=x+step;
                        weight=weight*exp(-SigmaMaj*OverSample*step*(1-s));
                        r=r-step*SigmaMaj*OverSample*s;
                        i_cell=i_cell+1;
                    else
                        x=x+xs;
                        weight=weight*exp(-SigmaMaj*OverSample*xs*(1-s));
                        inTravel=0;
                    end
                end
                weight=weight/s;
                Sigma=SigmaX(i_cell);
                %      x=x+xs;
                if(x>=a)
                    Inside=0;
                end
                if(x<a && rand<q*Sigma/(SigmaMaj*OverSample))
                    VirtCollision=0;
                    weight=weight/q; %*Sigma/SigmaMaj; % a q-s r?sz nem j?, valami nem ok? a s?lykorrekci?val
                    Contribution=x*weight;
                    ResWET=ResWET+Contribution;
                    VarWET=VarWET+Contribution^2;
                    % tally weight distribution
                    i=ceil(x/a*NWdiag);
                    ExpWeight(i)=ExpWeight(i)+weight;
                    VarWeight(i)=VarWeight(i)+weight^2;
                    NSample_per_bin(i)=NSample_per_bin(i)+1;
                    %tally scores for histogram making at specific bins
%                     if(i==i_HalfSigma)
%                         Sc_H=[Sc_H weight];
%                     end
                    if(i==i_OneSigma)
                   %     Sc_1=[Sc_1 weight];
                   Res_One(i_o,i_qs)=Res_One(i_o,i_qs)+weight;
                   Var_One(i_o,i_qs)=Var_One(i_o,i_qs)+weight^2;
                   Ns_One(i_o,i_qs)=Ns_One(i_o,i_qs)+1;
                    end
%                     if(i==i_TwoSigma)
%                         Sc_2=[Sc_2 weight];
%                     end
                    if(i==i_FiveSigma)
                  %      Sc_5=[Sc_5 weight];
                         Res_Five(i_o,i_qs)=Res_Five(i_o,i_qs)+weight;
                         Var_Five(i_o,i_qs)=Var_Five(i_o,i_qs)+weight^2;
                         Ns_Five(i_o,i_qs)=Ns_Five(i_o,i_qs)+1;
                    end
                end
                weight=weight*(OverSample*SigmaMaj-Sigma)/(OverSample*SigmaMaj-Sigma*q);
            end
        end

ResWET=ResWET/NSamp;

fprintf('Ideal = %f WET= %f relvarWET= %f\n', ResIdeal,ResWET);

ResWet=0;
    end
end

VarWeight=sqrt(VarWeight./(ExpWeight).^2-1./NSample_per_bin);
ExpWeight=ExpWeight./NSample_per_bin;
VarWET=sqrt(VarWET/ResWET^2-1/NSamp);
ResWET=ResWET/NSamp;
ResIdeal;
fprintf('Ideal = %f WET= %f relvarWET= %f\n', ResIdeal,ResWET,VarWET);

Var_One=sqrt((Var_One./(Res_One).^2-1./Ns_One).*Ns_One); % this is now the variance of the latent random variable not the variance of the mean of them! 
Var_Five=sqrt((Var_Five./(Res_Five).^2-1./Ns_Five).*Ns_Five);
%plot vs q
figure();
for i_o=1:N_o

    title('qZ=' +string(qZ)+ ' variances at One path length');
        xlabel('q');
        ylabel('relative variance');
        plot(q_a(i_o,:),Var_One(i_o, :));
        hold on;
end
legend(string(1:N_o));
hold off;
%plot vs s
figure();
for i_o=1:N_o
    
    title('qZ=' +string(qZ)+ 'variances at One path length');
        xlabel('s');
        ylabel('relative variance');
        plot(s_a(i_o,:),Var_One(i_o, :));
        hold on;
end
legend(string(1:N_o));
hold off;
%5-----------------
figure();
for i_o=1:N_o

    title('qZ=' +string(qZ)+ 'variances at Five path length');
        xlabel('q');
        ylabel('relative variance');
        plot(q_a(i_o,:),Var_Five(i_o, :));
        hold on;
end
legend(string(1:N_o));
hold off;
%plot vs s
figure();
for i_o=1:N_o
    
    title('qZ=' +string(qZ)+ 'variances at Five path length');
        xlabel('s');
        ylabel('relative variance');
        plot(s_a(i_o,:),Var_Five(i_o, :));
        hold on;
end
legend(string(1:N_o));
hold off;

% plotting
figure;
xx=1:NWdiag;
plot(xx,ExpWeight);
hold on;
%plot(xx,exp(-Sigma*xx*a/NWdiag)/NWdiag*a);
plot(xx,exp(-Sigma*xx/100*10*(1-qZ))/qZ);
%exp(-Sigma*OverSample*xs*(1-qZ))/qZ
title('Weight expectation - the ideal weight');
hold off;

figure;
xx=1:NWdiag;
plot(xx,VarWeight);
%hold on;
%plot(xx,exp(-Sigma*xx*(1-qZ))/qZ);
title('Relative Variance of Weights');
hold off;


% figure;
% histogram(Sc_H);
% title('At half mean free path');
% hold off;
% figure;
% histogram(Sc_1);
% title('At One mean free path');
% hold off;
% figure;
% histogram(Sc_2);
% title('At Two mean free paths');
% hold off;
% figure;
% histogram(Sc_5);
% title('At Five mean free paths');
% hold off;

