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
NS=5;
%NS=2;
%NS=5;
IF_O=2*IF_O/length(IF_O);
% HADTFD BASED
win_length=45;
%for snr=-10:2:10
iiii=0;
delta=2;
L=32*1;
FFT_length=length(Sig1);
snr=5;
NS=25;
mseee=ones(NS,3,128)*1;

for k1=1:NS
    a=2*rand(1,3);
    % a(:)=1;
    %SigO=(1+a(1))*Sig1+(1+a(2))*Sig2+(1+a(3))*Sig3;
    SigO=4*Sig1+2*Sig2+1*Sig3;
    
    %    SigO=Sig1+Sig2+Sig3;
    
    Sig=awgn(SigO,snr,'measured');
    
    findex =ADTFD_RANSAC(Sig,3,15,64,num,500/1,4,8,1);
    
    IF=zeros(1,length(Sig));
    dis=0;
    clear c;
    
    for ii=1:num
        
        t=1:SampFreq;
        IF=findex(ii,:)/length(Sig);
        %t=t(5:end-5);
        for i=1:num
            c(i)=sum(abs(IF(t)'-IF_O(t,i)).^2);
        end
        [a1, b1]=min(c);
        mseee(k1,b1,:)=abs(IF(t)'-IF_O(t,b1)).^2;
        if dis==1
            figure;
            plot(t,IF(t),'-',t,IF_O(t,b1),'d');
        end
    end
    
    
end



mse_ADTFD=mean(mseee);


t = 0:1/SampFreq:1-1/SampFreq;

plot(t, 10*(log10(reshape(mse_ADTFD(1,1,:),1,128))),'-rh','linewidth',4);
hold on;
plot(t, 10*(log10(reshape(mse_ADTFD(1,2,:),1,128))),'-bh','linewidth',4);
hold on;
plot(t, 10*(log10(reshape(mse_ADTFD(1,3,:),1,128))),'-kh','linewidth',4);
xlabel('Time (s)');
ylabel('Point wise mean equare error (dB)');

%legend('Strongest Component','2nd Strongest Component','Weakest Component');
legend('Component-1','Component-2','Component-3');

