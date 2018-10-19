n1=50; n2=10; n3=50; k=50; SNR=30;
A=rand(n1,n2,n3);

[Uc Sc Vc]=c_svd(A);


Ucr=Uc(:,1:n2,:);  %
Uc_r=Uc(:,n2+1:n1,:);  %

P_fa=zeros(SNR+1,1);

Pfa=0.2;



c=rand(n2,1,n3);
        yc0=c_product(Ucr,c);
        
        noise=randn(n1,1,n3);
        
        Yc0=0;
        
     for j=1:n3       
          ac=norm(yc0(:,1,j));
          bc=ac^2;
          Yc0=Yc0+bc;
     end
     clear ac bc;
     
for snr=1:SNR
    num=0;
    counter=0;
    
    SN=snr/10;
         
     sn=10^SN;
     
     yc=sqrt(n1*n3*sn)*yc0/sqrt(Yc0);
     ycn=yc+noise;
     
      while counter<20 && num<200 
           W(1:k)=randperm(n1,k);
        w=sort(W(1:k));
        
       Ycw=zeros(k,1,n3);
        Ucw=zeros(k,n2,n3);
        
        for i=1:k 
           Ycw(i,1,:)=ycn(w(i),1,:);
           Ucw(i,:,:)=A(w(i),:,:);  
        end
        
     %   Pw=zeros(k,k,n3);
        
      % Yw_c=my_dct(Ycw);
        %  Yw_b=bl_dia(Yw_c);
        
        %  Uw_c=my_dct(Ucw);
        %  Uw_b=bl_dia(Uw_c);
       %   a=Uw_b';
       a=c_transpos(Ucw);
       b=c_product(a,Ucw);
       f=my_dct(b);
        g=zeros(n2,n2,n3);
        for i=1:n3
        g(:,:,i)=pinv(f(:,:,i));
        end
        h=my_idct(g);
        s=c_product(Ucw,h);
        Pw=c_product(s,a);
        clear a b f g h s;
        
        a=c_product(Pw,Ycw);
        ec=Ycw-a;
        clear a ;
        
        n=0;
        for i=1:n3
            a=norm(ec(:,1,i));
            b=a^2;
            n=n+b;
        end
        
     mat_Pw=c_mat(Pw);
     r=rank(eye(k*n3)-mat_Pw);
     
     freedom=k*n3-r;
threshold=chi2inv(1-Pfa,freedom);
clear mat_Pw;
        
        if n>threshold
            counter=counter+1;
        end
        counter
        num=num+1
    end
    P_fa(snr+1)=counter/num;
end
  
snr=0:30;
figure;
plot(snr,P_fa(snr+1),'-o'); xlabel('snr'),ylabel('P_{fa}')
     

