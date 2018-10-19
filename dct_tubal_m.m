
clear;

n1=50; n2=10; n3=50; K=30; SNR=0;
A=rand(n1,n2,n3);

[Uc Sc Vc]=c_svd(A);


Ucr=Uc(:,1:n2,:);  %
Uc_r=Uc(:,n2+1:n1,:);  %


Pdc=zeros(K,1);w

Pfa=0.01;

c=rand(n1-n2,1,n3);
yc0=c_product(Uc_r,c);
        
noise=randn(n1,1,n3);



for k=11:K
    
   freedom=(k-n2)*n3;
   threshold=chi2inv(1-Pfa,freedom);
    
     num=0;
    counter=0;
    
    while counter<200&&num<2000
        
        Yc0=0;
        
     for j=1:n3       
          ac=norm(yc0(:,1,j));
          bc=ac^2;
          Yc0=Yc0+bc;
     end
     clear ac bc;
     
      yc=sqrt(n1*n3)*yc0/sqrt(Yc0);
      ycn=yc+noise;
        
        W(1:k)=randperm(n1,k);
        w=sort(W(1:k));
        
       Ycw=zeros(k,1,n3);
       Ucw=zeros(k,n2,n3);
        
        for i=1:k 
           Ycw(i,1,:)=ycn(w(i),1,:);
           Ucw(i,:,:)=A(w(i),:,:);  
        end
        
        Y=unfold(Ycw);
        U=c_mat(Ucw);
        a=U';
        b=a*U;
        d=pinv(b);
        f=U*d;
        Pw=f*a;
        
        clear a b d f;
        
        e=Y-Pw*Y;
        
        n=norm(e);
        n=n^2;
        
        if n<=threshold
            counter=counter+1;
        end
        counter
        num=num+1
    end
    Pdc(k)=(num-counter)/num;
    
end

k=11:K;
figure;
plot(k/n1,Pdc(k),'-o'); xlabel('Sampling rate'),ylabel('Pd'),title('Detection probability');
hold on;


save('data\dct_tubal_m_Pdc','Pdc');