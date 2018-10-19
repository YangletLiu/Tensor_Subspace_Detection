
clear;

n1=50; n2=10; n3=50; k=11;J=4;SNR=15;
A=rand(n1,n2,n3);
%A=round(rand(n1,n2,n3)*255);

[Uf Sf Vf]=t_svd(A);
%[Uc Sc Vc]=c_svd(A);

Ufr=Uf(:,1:n2,:);
Uf_r=Uf(:,n2+1:n1,:);

%Ucr=Uc(:,1:n2,:);  %
%Uc_r=Uc(:,n2+1:n1,:);  %


Pd=zeros(SNR+1,J);
P_fa=[0.1;0.01;0.001;0.0001];

figure;

for j1=1:J
    
    Pfa=P_fa(j1);

freedom=(k-n2)*n3;
threshold=chi2inv(1-Pfa,freedom);


for snr=0:SNR
    
    num=0;
    counter=0;
    
    SN=snr/10;
         
     sn=10^SN;
     
    
    while counter<500 && num<10000  
        
        c=rand(n1-n2,1,n3);
        yf0=t_product(Uf_r,c);
        noise=randn(n1,1,n3);
         Y0=0;
     for j=1:n3
         a=norm(yf0(:,1,j));
         b=a^2;
         Y0=Y0+b;
     end
     clear a b;
     
      yf=sqrt(n1*n3*sn)*yf0/sqrt(Y0);
        
     yfn=yf+noise;
        
        W(1:k)=randperm(n1,k);
        w=sort(W(1:k));
        Yw=zeros(k,1,n3);
        Uw=zeros(k,n2,n3);
        for i=1:k
        Yw(i,1,:)=yfn(w(i),1,:);
        Uw(i,:,:)=A(w(i),:,:);
        end
        a=transpos(Uw);
        b=t_product(a,Uw);
        f=fft(b,[],3);
        g=zeros(n2,n2,n3);
        for i=1:n3
        g(:,:,i)=pinv(f(:,:,i));
        end
        h=ifft(g,[],3);
        s=t_product(Uw,h);
        Pw=t_product(s,a);
        clear a b f g h s;
        a=t_product(Pw,Yw);
        e=Yw-a;
        n=0;
        for i=1:n3
            a=norm(e(:,1,i));
            b=a^2;
            n=n+b;
        end
        
        if n<=threshold
            counter=counter+1;
        end
        counter
        num=num+1
    end
    Pd(snr+1,j1)=(num-counter)/num;
    
end

snr=0:SNR;

plot(snr,Pd(snr+1,j1),'-o'); xlabel('snr'),ylabel('Pd'),title('Detection probability');
hold on;

end
legend('P_{fa}=10^{-1}','P_{fa}=10^{-2}','P_{fa}=10^{-3}','P_{fa}=10^{-4}');

save('data\fft_tubalt_niosy_k11_Pd','Pd');
