clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;

t = 0:1/SampFreq:1-1/SampFreq;

Sig1 = 1*exp(1i*(2*pi*(70*(t-0.5).^3))+1i*(2*pi*(0*t))); %300t»òÕß150t
Sig2 = 1*exp(1i*(2*pi*(70*(t-0.5).^3))+1i*(2*pi*(11*t))); %300t»òÕß150t

IF_O(:,1)=3*70*(t-0.5).^2/1;
IF_O(:,2)=3*70*(t-0.5).^2/1+11;


num=2;
NS=100;
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
iiii=iiii+1;

a=2*rand(1,3);
% a(:)=1;
%SigO=(1+a(1))*Sig1+(1+a(2))*Sig2+(1+a(3))*Sig3;
SigO=1*Sig1+1*Sig2;
delta=3;
%    SigO=Sig1+Sig2+Sig3;

Sig=awgn(SigO,10,'measured');
plot(t,IF_O/2,'g:','linewidth',3);
for kkkkk=0:2
hold on;    
    % ORIGINAL
    
    if kkkkk==0
        findex1 =FASTEST_IF(Sig,win_length, num, delta,L,0,0,2,length(Sig))*2*SampFreq;
        plot(t,findex1.'/256,'k--','linewidth',3);

    elseif kkkkk==1
        findex2 =ADTFD_RANSAC_Illust(Sig,3,15,64,num,500/1,4,8,1);
        plot(t,findex2.'/256,'b-.','linewidth',3);
    elseif kkkkk==2
        
        findex3 =QML_RANSAC(Sig,3,15,64,num,500/1,4,8,1);
        plot(t,findex3.'/256,'r:','linewidth',3);
        
    end
    
end
    %legend('Original IF','FAST-IF','ADTFD-RANSAC','QML-RANSAC');
    close all;
    figure;
    plot(t,IF_O/2,'g:','linewidth',3);

            findex2 =ADTFD_RANSAC_Illust(Sig,3,15,64,num,500/1,4,8,1);

    xlabel('Time(s)');
    ylabel('Instantaneous Frequency (Hz)');
I=HTFD_new1(Sig,3,8,64);
t = 0:1/SampFreq:1-1/SampFreq;
f=0:1/256:0.5-1/256;    
figure;
imagesc(t,f,I)
set(gcf,'Position',[20 100 640 500]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
%title('(c)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);



