n1=50; n2=10; n3=50; k=11;
A=rand(n1,n2,n3);

[Uf Sf Vf]=t_svd(A);
%[Uc Sc Vc]=c_svd(A);

Ufr=Uf(:,1:n2,:);
Uf_r=Uf(:,n2+1:n1,:);

%Ucr=Uc(:,1:n2,:);  %
%Uc_r=Uc(:,n2+1:n1,:);  %


P_fa=zeros(21,1);
Pfa=0.1;

freedom=(k-n2)*n3;
threshold=chi2inv(1-Pfa,freedom);

c=rand(n2,1,n3);
        yf0=t_product(Ufr,c);
        noise=randn(n1,1,n3);
         Y0=0;
     for j=1:n3
         a=norm(yf0(:,1,j));
         b=a^2;
         Y0=Y0+b;
     end
     clear a b;

for snr=0:10
    
    num=0;
    counter=0;
    
    SN=snr/10;
         
     sn=10^SN;
     
     yf=sqrt(n1*n3*sn)*yf0/sqrt(Y0);
        
     yfn=yf+noise;
    
    while counter<50 && num<1000  
        
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
        
        if n<threshold
            counter=counter+1;
        end
        counter
        num=num+1
    end
    P_fa(snr+1)=(num-counter)/num;
    
end

snr=0:9;
figure;
plot(snr,P_fa(snr+1),'-o'); xlabel('snr'),ylabel('P_{fa}'),title('Detection probability');
