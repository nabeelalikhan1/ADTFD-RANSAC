function F=ADTFD_RANSAC(Sig,a,b,WL,num,NS,T,number_samples,step)
for iii=1:num
[I,O]=HTFD_new1(Sig,a,b,WL);
%[I,O]= post_processing_directional(I,2,15,48);
%[I,O]= post_processing_directional(I,2,15,48);

O=O*3-90;
%figure; imagesc(I)
    [v,ind_a]=max(I);
    for ii=1:length(I)
    SI(:,ii)=sort(I(:,ii),'descend');
    end
    
    for i1=1:NS
        rr=randperm(length(Sig),length(Sig)/number_samples);
        IF=zeros(1,length(Sig));
        %for i=1:length(rr)
         %   IF(rr(i))=ind_a(rr(i));
        %end
            IF(rr(:))=ind_a(rr(:));

        %IF(1)=ind_a(1);
        %IF(end)=ind_a(end);
        %tfsapl(Sig,I);
        %imagesc(I)
        ind=find(IF>0);
       % IF=IF_param(IF(:),1,ind);
                IF=IF_param(IF(:),1,ind);

        %IF=interp1(ind,IF(ind),1:length(Sig),'linear','extrap');
        %figure; plot(IF)
        IF=round(IF);
        IF(IF<1)=1;
        IF(IF>length(Sig))=length(Sig);
        for i=1:length(IF)
           if and(IF(i)-T>0,IF(i)+T<length(Sig))
               [~,aa]=max(I(IF(i)-T:IF(i)+T,i));
               IF(i)=IF(i)-T+aa-1;
           elseif IF(i)+T>=length(Sig)
               [~,aa]=max(I(IF(i)-T:end,i));
               IF(i)=IF(i)-T+aa-1;
           else
               [~,aa]=max(I(1:T,i));
               IF(i)=aa;
               
           end
            
                        v=find(SI(:,i)==I(IF(i),i));

            
            IA(i)=v(1)-1;
            IO(i)=O(IF(i),i);
            
           % II(IF(i),i)=I(IF(i),i);
        end
       % IF=IF_param(IF(:),2,1:length(I));
            IFF(i1,:)=IF;
     %m(i1)=1*sum(IA)+1*sum(abs(diff(IO(2:8:end))))+1*sum(abs(diff(IF(1:8:end))))+0*sum(abs(diff(IF,2)));
            m(i1)=1*sum(IA(1:step:end))+1*sum(abs(diff(IO(1:step:end))))+1*sum(abs(diff(IF(1:step:end))))+0*sum(abs(diff(IF(1:step:end),2)));
%             IF=IF/(2*length(Sig));
%         Phase=2*pi*filter(1,[1 -1],IF);
%         s_dechirp=exp(-1i*Phase).';
%         m(i1)=sum(abs(Sig.*s_dechirp));
    end
        
        [~,indd]=min(m);
        F(iii,:)=IFF(indd,:);
        IF=F(iii,:);
        IF=IF/(2*length(Sig));
        Phase=2*pi*filter(1,[1 -1],IF);
        s_dechirp=exp(-1i*Phase);
        
        
        LL=T;
        %TF filtering for each sensor
        s1 = Sig.*(s_dechirp);
        s2=fftshift(fft(s1));
        %figure; plot(abs(s2));
            %  e_max
           % if sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2))>e_max
            %    e_max=sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2));
            %end
            % Energy of the last component
            s2(length(Sig)/2-LL:length(Sig)/2+LL)=0;
            s2=ifft(ifftshift(s2)).*conj(s_dechirp);
            Sig=s2;%-extr_Sig(iii);
        
        
end
end
       