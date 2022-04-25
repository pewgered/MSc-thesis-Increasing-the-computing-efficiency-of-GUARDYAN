% IMPORTANT!
% CHANGE OF NOTATIIONS: NOW SIGMA SAMP IS THE TOTAL SAMPLING XS, THAT IS IT
% IS THE PRODUCT OF SIGMA MAJ, OVERSAMPLING AND THE S STRETCHING PARAMETER.
% THE IDEA IS TO LATER DECIDE WHICH PART OF SIGMA S IS USED FOR STRETCHING
% AND WHICH FOR COLLISION SAMPLING
clc
clear all
close all
% integral of versions of x exp(-x) as I lost the original
% Woodcock exponential transform test
% this version tries to set the Woodcock miss-weight to zero
N_conv_step=10;
NSamp=1E4;
a=10;
Diff=0.2;
ResWET(1:N_conv_step)=0;
VarWET(1:N_conv_step)=0;
NWdiag=100;
N_Adj_Res=20;
N_qs=1;
N_o=1;
ExpWeight(1:N_conv_step,1:NWdiag)=0;
VarWeight(1:N_conv_step,1:NWdiag)=0;
NSample_per_bin(1:N_conv_step,1:NWdiag)=0;
for i_conv=1:N_conv_step
    N_Adj_Res=i_conv*10;
    Sigma=1.0;
    SigmaX(1:N_Adj_Res+1)=1.0;
    SigmaZ(1:N_Adj_Res+1)=0;
    SigmaS(1:N_Adj_Res+1)=0;
    for i=1:N_Adj_Res+1
        x_coo=a/N_Adj_Res*(i-1)+a/N_Adj_Res/2; %to calculate the adjoint at interval middles
        SigmaZ(i)=x_coo*exp(-x_coo)/(exp(-x_coo)*(x_coo+1)-exp(-a)*(a+1));
        
        %    if (mod(i,2)==0)
        SigmaS(i)=max(Sigma,SigmaZ(i))*1.0; % the default oversampling is to always have virtual collisions
        %    if (mod(i,2)==0)
        %        SigmaS(i)=SigmaS(i)*2.0;
        %    end
    end
    ResIdeal=1-exp(-a)*(a+1);
    
    % expected value and distributiontest for the Woodcock Exponential
    % Transform
    qZ=2.0; % the ideal multiplication foactor to Sigma_local Sigma_z=Sigma(x)*qz
    
    q_a(N_o,1:N_qs)=0;
    s_a(N_o,1:N_qs)=0;
    
    for i_o=1:N_o
    for i_qs=1:N_qs
    OverSample=i_o; % oversampling is set in a fashion that it does not affect the expected value of the cross section sampled at the end
    % for debug

    % virtual/real collision selsection chance
    % max possible q to stay below 1 is:
    SigmaMaj=SigmaS(1); % this is supposed to account for assuming a lack of a tight majorant at every phase space loation
    % special spots for histogram. Any further modification to this value
    % happens through the oversampling parameter. SigmaMaj=Sigma would make Exp
    % Transformation alone the ideal choice.
    %         q_max=SigmaMaj*OverSample/Sigma;
    %         q_max=1.8; % now zooming in
    %         q_min=0.05; % just not to make it 0
    %         q_min=1.0; % now zooming in
    %         q=(q_max-q_min)/N_qs*i_qs+q_min;
    % ideal q is set now. For later simulations it must move to later!

    q=qZ*(1+Sigma/(SigmaMaj*OverSample)*(1-qZ)  );
    q_a(i_o,i_qs)=q;
    %q=qZ;
    % %for debug
    % s=qZ;
    % q=Sigma/(SigmaMaj*OverSample/s);
    %
    % %
    s=qZ/q; % this ensures that the ideal qz is sampled

    if(q*Sigma/(SigmaMaj*OverSample)>1)
        q=1.0 % this is left here for warning
    end

    s_a(i_o,i_qs)=s;
    % exponential transform parameter
    %     fprintf('qz = %f o= %f q= %f s= %f\n', qZ,OverSample,q,s);

    for i=1:NSamp
        VirtCollision=1;
        weight=1.0;
        Inside=1;
        x=0;
        i_cell=1;
        while(VirtCollision && Inside)
            rr=rand;
            r=-log(1-rr); %optical distance
            inTravel=1;
            while(inTravel && i_cell<=N_Adj_Res)
                SigmaMaj=SigmaS(i_cell);
                % to set the ideal paramters
                qZ=SigmaZ(i_cell)/Sigma;
                % q=qZ*(1+Sigma/(SigmaMaj*OverSample)*(1-qZ)  ); % Now this is the new stuff!
                %q_a(i_o,i_qs)=q;
                q=1;
                s=qZ/q; % this ensures that the ideal qz is sampled

                %s should also depend on i_cell later
                %xs=r/(SigmaMaj*OverSample*s); %
                xs=r/SigmaS(i_cell);
                
                if(x+xs>=a/N_Adj_Res*i_cell)
                    step=a/N_Adj_Res*i_cell-x;
                    x=x+step;
                    %weight=weight*exp(-SigmaMaj*OverSample*step*(1-s)); %
                    %r=r-step*SigmaMaj*OverSample*s;%
                    r=r-step*SigmaS(i_cell);
                    i_cell=i_cell+1;
                else
                    x=x+xs;
                    %weight=weight*exp(-SigmaMaj*OverSample*xs*(1-s));%
                    weight=exp(-x)/(SigmaS(i_cell)*(1-rr));
                    inTravel=0;
                end
            end
            %weight=weight/s;
            Sigma=SigmaX(i_cell);
            %      x=x+xs;
            if(x>=a)
                Inside=0;
            end
            if(x<a && rand<q*Sigma/(SigmaMaj*OverSample)) %
               VirtCollision=0;
                weight=weight/q; %
                Contribution=x*weight;
                ResWET(i_conv)=ResWET(i_conv)+Contribution;
                VarWET(i_conv)=VarWET(i_conv)+Contribution^2;
                % tally weight distribution
                i=ceil(x/a*NWdiag);
                ExpWeight(i_conv,i)=ExpWeight(i_conv,i)+weight;
                VarWeight(i_conv,i)=VarWeight(i_conv,i)+weight^2;
                NSample_per_bin(i_conv,i)=NSample_per_bin(i_conv,i)+1;
                %tally scores for histogram making at specific bins
                %                     if(i==i_HalfSigma)
                %                         Sc_H=[Sc_H weight];
                %                     end
            end
            weight=weight*(OverSample*SigmaMaj-Sigma)/(OverSample*SigmaMaj-Sigma*q); %
        end
    end
    end
    end
    VarWeight(i_conv,:)=sqrt(VarWeight(i_conv,:)./(ExpWeight(i_conv,:)).^2-1./NSample_per_bin(i_conv,:));
    ExpWeight(i_conv,:)=ExpWeight(i_conv,:)./NSample_per_bin(i_conv,:);
    VarWET(i_conv)=sqrt(VarWET(i_conv)/ResWET(i_conv)^2-1/NSamp);
    ResWET(i_conv)=ResWET(i_conv)/NSamp;
    ResIdeal;
    fprintf('Ideal= %f WET= %f relvarWET= %f NumOfIntervals= %d\n', ResIdeal,ResWET(i_conv),VarWET(i_conv), 10*i_conv);
end

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