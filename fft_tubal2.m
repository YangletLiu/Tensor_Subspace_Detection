n1=50; n2=10; n3=50; T=100;
A=rand(n1,n2,n3);
[U S V]=t_svd(A);
Ur=U(:,1:n2,:);
c=rand(n2,1,n3);
X=t_product(Ur,c);
clear c;
X_c=fft(X,[],3);

X_b=bl_dia(X_c);


X0=norm(X_b, 'fro');

a=transpos(Ur);
Ps=t_product(Ur,a);
clear a;
U_r=U(:,n2+1:n1,:);
c=rand(n1-n2,1,n3);
y=t_product(U_r,c);
clear c;


Y_c=fft(y,[],3);

Y_b=bl_dia(Y_c);


Y0=norm(Y_b,'fro');



E_Y0=zeros(n1,T);
Es_X0=zeros(n1,T);    

for t=1:T
W=zeros(n1,1);

E=zeros(n1,1);
Es=zeros(n1,1);

for k=1:n1
W(1:k)=randperm(n1,k);
w=sort(W(1:k));
Yw=zeros(k,1,n3);
Xw=zeros(k,1,n3);
Uw=zeros(k,n2,n3);
for i=1:k
Yw(i,1,:)=y(w(i),1,:);
Xw(i,1,:)=X(w(i),1,:);
Uw(i,:,:)=A(w(i),:,:);

end

Yw_f=fft(Yw,[],3);
Yw_b=bl_dia(Yw_f);
Xw_f=fft(Xw,[],3);
Xw_b=bl_dia(Xw_f);

Uw_f=fft(Uw,[],3);
Uw_b=bl_dia(Uw_f);
a=Uw_b';
b=a*Uw_b;
c=pinv(b);
d=Uw_b*c;
Pw_b=d*a;
clear a b c d w Xw Yw Uw Uw_c Uw_b;
e=Yw_b-Pw_b*Yw_b;
es=Xw_b-Pw_b*Xw_b;
     
E(k)=norm(e,'fro');
Es(k)=norm(es,'fro');
      
E(k)

end

E_Y0(:,t)=E./Y0;
Es_X0(:,t)=Es./X0;    
        %
end
Lmean=zeros(n1,1);
Lmax=zeros(n1,1);
Lmin=zeros(n1,1);
Lsmean=zeros(n1,1);
Lsmax=zeros(n1,1);
Lsmin=zeros(n1,1);

for i=1:n1

 Lmean(i)=mean(E_Y0(i,:));
Lmax(i)=max(E_Y0(i,:));
Lmin(i)=min(E_Y0(i,:));
Lsmean(i)=mean(Es_X0(i,:));
Lsmax(i)=max(Es_X0(i,:));
Lsmin(i)=min(Es_X0(i,:));

end

N=zeros(n1,1);
for j=1:n1
n=0;
for i=1:n1
a=0;
for k=1:n3
h=Ps(i,j,k)^2;
a=a+h;
end
n=n+a;
end
N(j)=n;
end
c=max(N);
mu=n1*c/n2;
mu

k=1:n1;
figure;
plot(k,Lmin,'x',k,Lmax,'+',k,Lmean,'*');axis([0 n1 0 1.2]);xlabel('m'),title('Projection Residual');
legend('minimum','maximum','mean');
grid minor;
figure;
plot(k,Lsmin,'x',k,Lsmax,'+',k,Lsmean,'*');axis([0 n1 0 0.2]);xlabel('m'),title('Projection Residual');
legend('minimum','maximum','mean');
grid minor;