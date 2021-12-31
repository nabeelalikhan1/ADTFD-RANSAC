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
for snr=0:2:10%snr=-10:2:10
    iiii=iiii+1;
    
    for k1=1:NS
        a=2*rand(1,3);
       % a(:)=1;
        %SigO=(1+a(1))*Sig1+(1+a(2))*Sig2+(1+a(3))*Sig3;
            %SigO=(1+a(1))*Sig1+(1+a(2))*Sig2;

        SigO=4*Sig1+2*Sig2+1*Sig3;

        Sig=awgn(SigO,snr,'measured');
        
        for kkkkk=0:11
            
            % ORIGINAL
            switch kkkkk
                case 0
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,50/1,4,8,1);
                case 1
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,100/1,4,8,1);
                case 2
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,200/1,4,8,1);
                case 3
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,500/1,4,8,1);
                case 4
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,50/1,4,4,1);
                case 5
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,100/1,4,4,1);
                case 6
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,200/1,4,4,1);
                case 7
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,500/1,4,4,1);
                case 8
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,50/1,4,2,1);
                case 9
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,100/1,4,2,1);
                case 10
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,200/1,4,2,1);
                case 11
                    findex =ADTFD_RANSAC(Sig,3,15,64,num,500/1,4,2,1);
                    
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
            
            switch kkkkk
                
                case 0
                   mse_250_12(k1)=mean(msee);
                case 1
                   mse_500_12(k1)=mean(msee);
                case 2
                   mse_1000_12(k1)=mean(msee);
                case 3
                   mse_2000_12(k1)=mean(msee);
                case 4
                   mse_250_25(k1)=mean(msee);
                case 5
                   mse_500_25(k1)=mean(msee);
                case 6
                  mse_1000_25(k1)=mean(msee);

                case 7
                  mse_2000_25(k1)=mean(msee);

                case 8
                   mse_250_50(k1)=mean(msee);
                case 9
                   mse_500_50(k1)=mean(msee);
                case 10
                  mse_1000_50(k1)=mean(msee);

                case 11
                  mse_2000_50(k1)=mean(msee);
            
        end
        end
    end
                   mmse_250_12(iiii)=mean(mse_250_12);
                   mmse_500_12(iiii)=mean(mse_500_12);
                   mmse_1000_12(iiii)=mean(mse_1000_12);
                   mmse_2000_12(iiii)=mean(mse_2000_12);
                   mmse_250_25(iiii)=mean(mse_250_25);
                   mmse_500_25(iiii)=mean(mse_500_25);
                   mmse_1000_25(iiii)=mean(mse_1000_25);
                   mmse_2000_25(iiii)=mean(mse_2000_25);
                   mmse_250_50(iiii)=mean(mse_250_50);
                   mmse_500_50(iiii)=mean(mse_500_50);
                   mmse_1000_50(iiii)=mean(mse_1000_50);
                   mmse_2000_50(iiii)=mean(mse_2000_50);
        
        
    
    end
    
    figure;
    snr=0:2:10;
    plot(snr, 10*(log10(mmse_250_12)),'-.ro','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_500_12)),'-.bo','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_1000_12)),'-.go','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_2000_12)),'-.ko','linewidth',4);
       hold on;
    plot(snr, 10*(log10(mmse_250_25)),'r:+','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_500_25)),'g:+','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_1000_25)),'b:+','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_2000_25)),'k:+','linewidth',4);
       hold on;
    plot(snr, 10*(log10(mmse_250_50)),'r--*','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_500_50)),'g--*','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_1000_50)),'b--*','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_2000_50)),'k--*','linewidth',4);
   
    xlabel('Signal to Noise Ratio');
    ylabel('Mean Square Error (dB)');
    legend('K=50,No=16','K=100,No=16','K=200,No=16','K=500,No=16','K=50,No=32','K=100,No=32','K=200,No=32','K=500,No=32','K=50,No=64','K=100,No=64','K=200,No=64','K=500,No=64');


    
figure;
    plot(snr, 10*(log10(mmse_250_12)),'ro:','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_500_12)),'bo:','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_1000_12)),'go:','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_2000_12)),'ko:','linewidth',4);
       hold on;
    plot(snr, 10*(log10(mmse_250_25)),'r+-.','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_500_25)),'g+-.','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_1000_25)),'b+-.','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_2000_25)),'k+-.','linewidth',4);
       hold on;
    plot(snr, 10*(log10(mmse_250_50)),'rs--','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_500_50)),'gs--','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_1000_50)),'bs--','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mmse_2000_50)),'ks--','linewidth',4);
   
    xlabel('Signal to Noise Ratio');
    ylabel('Mean Square Error (dB)');
    legend('K=50,No=16','K=100,No=16','K=200,No=16','K=500,No=16','K=50,No=32','K=100,No=32','K=200,No=32','K=500,No=32','K=50,No=64','K=100,No=64','K=200,No=64','K=500,No=64');

    