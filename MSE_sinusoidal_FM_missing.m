clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;
Sig1=exp(1*1i*(2*pi*(SampFreq*t/4 +2*sin(2*pi*t))));
Sig2=exp(1*1i*(2*pi*(SampFreq*t/4 -2*sin(2*pi*t))));
%Sig=exp(1*1i*(2*pi*(SampFreq*t/4 +2*sin(2*pi*t))))+exp(1*1i*(2*pi*(SampFreq*t/4 -2*sin(2*pi*t))))+1*0.*exp(1i*2*pi*SampFreq*t/4);
IF_O(:,1)=2*cos(2*pi*t)*2*pi+SampFreq/4;
IF_O(:,2)=-2*cos(2*pi*t)*2*pi+SampFreq/4;


%Sig=Sig.*([1:128 128:-1:1]);
num=2;
NS=250/25;
IF_O=2*IF_O/length(IF_O);
% HADTFD BASED
win_length=45;
%for snr=-10:2:10
iiii=0;
delta=2;
L=32*2;
FFT_length=length(Sig1);
for snr=16:16:80
    iiii=iiii+1;
    
    for k1=1:NS
        a=rand(1,2);
       % a(:)=1;
        Sig=(1+a(1))*Sig1+(1+a(2))*Sig2;
        Sig(randperm(length(Sig),snr))=0;
        
        for kkkkk=0:2
            
            % ORIGINAL
            
            if kkkkk==0
                findex =FASTEST_IF(Sig,win_length, num, delta,L,0,0,2,length(Sig))*2*SampFreq;
                
            elseif kkkkk==1
%                                        findex =ADTFD_RANSAC(Sig,3,15,64,num,500/1,4,32,8);
                                      %  findex =ADTFD_RANSAC(Sig,3,15,64,num,500/1,4,16/2,1);
                                                 findex =ADTFD_RANSAC_new(Sig,3,15,64,num,500/1,4,16/2,1,1);

            elseif kkkkk==2
                                       % findex =ADTFD_RANSAC1(Sig,3,15,64,num,500/1,4,16/2,1);
                                    %    findex =ADTFD_RANSAC1(Sig,3,15,64,num,500/1,4,16/2,2);
                                          findex =QML_RANSAC(Sig,3,15,64,num,500/1,4,16/2,1);

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

    snr=16:16:80;
    plot(snr, 10*(log10(mse_FAST_IF)),'-rh','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mse_R)),'-bh','linewidth',4);
    hold on;
    plot(snr, 10*(log10(mse_QML_RANSAC)),'-kh','linewidth',4);
    xlabel('Signal to Noise Ratio');
    ylabel('Mean Square Error (dB)');
    legend('FAST-IF','ADTFD-RANSAC','QML-RANSAC');
    

ylabel('Number of missing samples');
legend('FAST-IF','ADTFD-RANSAC','QML-RANSAC');
% hold on;
% plot(snr, 10*(log10(mse_proposed)),'-.k+','linewidth',4);
%
% hold on;
% plot(snr, 10*(log10(mse_mb)),'-.y+','linewidth',4);
%
% hold on;
% plot(snr, 10*(log10(mse_spec)),'-.g+','linewidth',4);
%



%xlabel('Signal to Noise Ratio');
%ylabel('Mean Square Error (dB)');
%legend('FAST_IF','Proposed_IF_16','Proposed_IF_8','Proposed_IF_4');
% axis([min(snr) max(snr)  -50  0])