clear;

n1=50; n2=10; n3=50; snr=0;J=4;K=30;
A=rand(n1,n2,n3);

[Uf Sf Vf]=t_svd(A);

Ufr=Uf(:,1:n2,:);
Uf_r=Uf(:,n2+1:n1,:);

Pd=zeros(K,J);
P_fa=[0.1;0.01;0.001;0.0001];

 %SN=snr/10;
 %sn=10^SN;
 figure;

 for j1=1:J
    
    Pfa=P_fa(j1);

for k=11:K
    
     num=0;
    counter=0;
    
    freedom=(k-n2)*n3;
    threshold=chi2inv(1-Pfa,freedom);
    
     while counter<100&&num<1000
         
          c=rand(n1-n2,1,n3);
        yf0=c_product(Uf_r,c);
        noise=randn(n1,1,n3);
        Yf0=0;
        for j=1:n3       
           ac=norm(yf0(:,1,j));
           bc=ac^2;
           Yf0=Yf0+bc;
        end
        clear ac bc; 
        yf=sqrt(n1*n3)*yf0/sqrt(Yf0);
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
    Pd(k,j1)=(num-counter)/num;
    
 end
    
 k=11:K;

plot(k*n3,Pd(k,j1),'-o'); xlabel('m'),ylabel('P_d'),title('Detection probability');
hold on;
end   
    
lengend('P_{fa}=10^{-1}','P_{fa}=10^{-2}','P_{fa}=10^{-3}','P_{fa}=10^{-4}');

save('data\fft_tubal_m_Pd','Pd');

%load E:\Li\tensor\matlab\msd_matlab\-mat\fft_tubal_snr_fa_Pd.mat;
%%导入将文件中的全部变量
