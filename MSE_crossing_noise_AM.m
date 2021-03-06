clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;

Sig1 = exp(-abs(t-0.5)/2).*exp(1i*(2*pi*(70*(t-0.5).^3))+1i*(2*pi*(0*t))); %300t»òÕß150t
Sig2 = exp(-abs(t-0.5)/2).*exp(1i*(-2*pi*(80*(t-0.5).^3))+1i*(2*pi*(64*t))); %300t»òÕß150t
Sig3 = exp(-abs(t-0.5)/2).*exp(1i*(2*pi*(70*(t-0.5).^3))+1i*(2*pi*(15*t))); %300t»òÕß150t

IF_O(:,1)=3*70*(t-0.5).^2/1;
IF_O(:,2)=-3*80*(t-0.5).^2/1+64;
IF_O(:,3)=3*70*(t-0.5).^2/1+15;


num=3;
NS=100;
%NS=2;
%NS=2;
%NS=5;
IF_O=2*IF_O/length(IF_O);
% HADTFD BASED
win_length=45;
%for snr=-10:2:10
iiii=0;
delta=3;
L=32*1;
FFT_length=length(Sig1);
for snr=0:2:20%snr=-10:2:10
    iiii=iiii+1;
    
    for k1=1:NS
        a=2*rand(1,3);
       % a(:)=1;
        %SigO=(1+a(1))*Sig1+(1+a(2))*Sig2+(1+a(3))*Sig3;
                SigO=4*Sig1+2*Sig2+1*Sig3;

            %    SigO=Sig1+Sig2+Sig3;

        Sig=awgn(SigO,snr,'measured');
        
        for kkkkk=0:2
            
            % ORIGINAL
            
            if kkkkk==0
                findex =FASTEST_IF(Sig,win_length, num, delta,L,0,0,2,length(Sig))*2*SampFreq;
                
            elseif kkkkk==1
                               findex =ADTFD_RANSAC(Sig,3,15,64,num,500/1,4,8,1);

            elseif kkkkk==2
                                       
 findex =QML_RANSAC(Sig,3,15,64,num,500/1,4,8,1);


            end  
            msee=0.1*ones(1,num);
            IF=zeros(1,length(Sig));
            dis=0;
            clear c;
            
            for ii=1:num
                
                t=1:SampFreq;
                IF=findex(ii,:)/length(Sig);
                t=t(5:end-5);
                for i=1:num
                    c(i)=sum(abs(IF(t)'-IF_O(t,i)).^2);
                end
                [a1 b1]=min(c);
                if msee(b1)>=a1(1)/length(Sig)
                    msee(b1)=a1(1)/length(Sig);
                end
                if dis==1
                    figure;
                    plot(t,IF(t),'-',t,IF_O(t,b1),'d');
                end
            end 
                if kkkkk==0
                    mse_FAST_IF_1(k1)=mean(msee);
                elseif kkkkk==1
                    mse_RANSAC(k1)=mean(msee);
                    elseif kkkkk==2
                    mse_RANSAC1(k1)=mean(msee);
                
            end
            
        end
        
        
     
        mse_FAST_IF(iiii)=mean(mse_FAST_IF_1);
        mse_R(iiii)=mean(mse_RANSAC);
        mse_QML_RANSAC(iiii)=mean(mse_RANSAC1);
        
    end
     mse_FAST_IF(iiii)
     mse_R(iiii)
     mse_QML_RANSAC(iiii)
    end  

    snr=0:2:20;
    plot(snr, 10*(log10(mse_FAST_IF)),'-rh','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mse_R)),'-bh','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mse_QML_RANSAC)),'-kh','linewidth',4);
    xlabel('Signal to Noise Ratio');
    ylabel('Mean Square Error (dB)');
    legend('FAST-IF','ADTFD-RANSAC','QML-RANSAC');
    
