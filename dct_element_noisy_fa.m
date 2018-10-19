
clear;

n1=50; n2=10; n3=50; SNR=15; %����
m=n1*n3;n=n2*n3;
k=11*n3;

A=rand(n1,n2,n3);

[Uc Sc Vc]=c_svd(A);


Ucr=Uc(:,1:n2,:);  %
Uc_r=Uc(:,n2+1:n1,:);  %

Pd=zeros(SNR+1,1);
Pfa=0.1;


freedom=k-n2*n3;
threshold=chi2inv(1-Pfa,freedom);

A_c=c_mat(A);
     
W=zeros(m,1);
     
for snr=0:SNR
    
    num=0;
    counter=0;
    
    SN=snr/10;
         
     sn=10^SN;
     
    
    while counter<50 && num<1000  
     %   
     c=rand(n2,1,n3);
        yc0=c_product(Ucr,c);
        noise=randn(n1,1,n3);
         Y0=0;
     for j=1:n3
         a=norm(yc0(:,1,j));
         b=a^2;
         Y0=Y0+b;
     end
     clear a b;
     
      yc=sqrt(n1*n3*sn)*yc0/sqrt(Y0);
        
     ycn=yc+noise;
     
     Y_u=unfold(ycn);
     
        W(1:k)=randperm(m,k);
        w=sort(W(1:k));
        Yw=zeros(k,1);
        Uw=zeros(k,n);
        
        for i=1:k
            Yw(i)=Y_u(w(i));
         
            Uw(i,:)=A_c(w(i),:);
           
        end
        %����Pw
        a=Uw';
        b=a*Uw;
        f=pinv(b);
        g=Uw*f;
        Pw=g*a;
        
        clear a b f g;
        
        a=Pw*Yw;
        
        e=Yw-a;
        
        clear a  ;
        %��������
        a=norm(e);
        E=a^2;
  %      
        if E>threshold
            counter=counter+1;
        end
        counter
        num=num+1
    end
    Pd(snr+1)=counter/num;
    
end

snr=0:SNR;
figure;
plot(snr,Pd(snr+1),'-o'); xlabel('snr'),ylabel('P_{fa}'),title('False alarm probability');