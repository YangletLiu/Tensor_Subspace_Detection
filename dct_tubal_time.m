n1=50; n2=10; n3=50; T=1;

A=rand(n1,n2,n3); 

[U S V]=c_svd(A);
Ur=U(:,1:n2,:);  %
U_r=U(:,n2+1:n1,:);  %

Ur_c=my_dct(Ur);
%U_r_c=my_dct(U_r);
Ur_b=bl_dia(Ur_c);
%U_r_b=bl_dia(U_r_c);
a=Ur_b';
b=a*Ur_b;
c=pinv(b);
d=Ur_b*c;
P_b=d*a;

a=c_transpos(Ur);

P=c_product(Ur,a);

N=zeros(n1,1);
for j=1:n1
n=0;
for k=1:n3
h=norm(P(:,j,k));
a=h^2;
n=n+a;
end
N(j)=n;
end
u=max(N);
mu=n1*u/n2;
clear N n a h ;

C=rand(n2,1,n3);
%C_c=my_dct(C);
%C_b=bl_dia(C_c);

X=c_product(Ur,C);
clear C;
C=rand(n1-n2,1,n3);
Y=c_product(U_r,C);
clear C;

X_c=my_dct(X);
Y_c=my_dct(Y);
X_b=bl_dia(X_c);
Y_b=bl_dia(Y_c);

X0=norm(X_b, 'fro');
Y0=norm(Y_b,'fro');

clear a b c d;

E_Y0=zeros(n1,T);
Es_X0=zeros(n1,T);    

for t=1:T
    W=zeros(n1,1);
    E=zeros(n1,1);
    Es=zeros(n1,1);   
   
   for k=15:15
       W(1:k)=randperm(n1,k);
       w=sort(W(1:k));
       Yw=zeros(k,1,n3);
       Xw=zeros(k,1,n3);
       Uw=zeros(k,n2,n3);
       for i=1:k
           Yw(i,1,:)=Y(w(i),1,:);
           Xw(i,1,:)=X(w(i),1,:);
           Uw(i,:,:)=A(w(i),:,:);
          
       end
      Yw_c=my_dct(Yw);
      Yw_b=bl_dia(Yw_c);
      Xw_c=my_dct(Xw);
      Xw_b=bl_dia(Xw_c);
     
      tic;
      Uw_c=my_dct(Uw);
      Uw_b=bl_dia(Uw_c);
      a=Uw_b';
      b=a*Uw_b;
      c=pinv(b);
      d=Uw_b*c;
      Pw_b=d*a;
      toc;
      clear a b c d w Xw Yw Uw Uw_c Uw_b;
      tic;
      e=Yw_b-Pw_b*Yw_b;
      es=Xw_b-Pw_b*Xw_b;
     
      E(k)=norm(e,'fro');
      Es(k)=norm(es,'fro');
      toc;
      E(k);
   end
   E_Y0(:,t)=E./Y0;
    Es_X0(:,t)=Es./X0;    
   
end

