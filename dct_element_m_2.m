
clear;

n1=50; n2=10; n3=50; snr=0;J=4; K=30;%参数
m=n1*n3;n=n2*n3;


A=rand(n1,n2,n3);

[Uc Sc Vc]=c_svd(A);


Ucr=Uc(:,1:n2,:);  %
Uc_r=Uc(:,n2+1:n1,:);  %

Pd=zeros(K,J);
P_fa=[0.1;0.01;0.001;0.0001];

figure;

 j1=2
    
    Pfa=P_fa(j1);



A_c=c_mat(A);
     
W=zeros(m,1);
     
for k=11:K
    
    freedom=(k-n2)*n3;
threshold=chi2inv(1-Pfa,freedom);

     num=0;
     counter=0;
     k_n3=k*n3; 
    
    while counter<100 && num<1000  
     %   
     c=rand(n1-n2,1,n3);
        yc0=c_product(Uc_r,c);
        noise=randn(n1,1,n3);
         Y0=0;
     for j=1:n3
         a=norm(yc0(:,1,j));
         b=a^2;
         Y0=Y0+b;
     end
     clear a b;
     
      yc=sqrt(n1*n3)*yc0/sqrt(Y0);
        
     ycn=yc+noise;
     
     Y_u=unfold(ycn);
     
        W(1:k_n3)=randperm(m,k_n3);
        w=sort(W(1:k_n3));
        Yw=zeros(k_n3,1);
        Uw=zeros(k_n3,n);
        
        for i=1:k_n3
            Yw(i)=Y_u(w(i));
         
            Uw(i,:)=A_c(w(i),:);
           
        end
        %计算Pw
        a=Uw';
        b=a*Uw;
        f=pinv(b);
        g=Uw*f;
        Pw=g*a;
        
        clear a b f g;
        
        a=Pw*Yw;
        
        e=Yw-a;
        
        clear a  ;
        %计算能量
        a=norm(e);
        E=a^2;
  %      
        if E<=threshold
            counter=counter+1;
        end
        counter
        num=num+1
    end
    Pd(k,j1)=(num-counter)/num;
    
end

k=0:K;

plot(k*n3,Pd(k,j1),'-o'); xlabel('m'),ylabel('P_d'),title('Detection probability');
hold on;



%  lengend('P_{fa}=10^{-1}','P_{fa}=10^{-2}','P_{fa}=10^{-3}','P_{fa}=10^{-4}');

save('data\dct_element_m_2_Pd','Pd');