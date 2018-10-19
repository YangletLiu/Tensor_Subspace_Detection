n1=50; n2=10; n3=50; T=100;
m=n1*n3; n= n2*n3;

A=rand(n1,n2,n3); 

[U S V]=c_svd(A);
Ur=U(:,1:n2,:);  %
U_r=U(:,n2+1:n1,:);  %

%U_c=c_mat(Ur);


a=c_transpos(Ur);

P=c_product(Ur,a);

N=zeros(n1,1);
for j=1:n1
n0=0;
for k=1:n3
h=norm(P(:,j,k));
a=h^2;
n0=n0+a;
end
N(j)=n0;
end
u=max(N);
mu=n1*u/n2;
clear N n0 a h ;

A_c=c_mat(A);

a=A_c';
b=a*A_c;
f=pinv(b);
g=A_c*f;
Ps=g*a;
clear a b f g;
N=zeros(m,1);
for j=1:n
    n0=norm(Ps(:,j));
    N(j)=n0^2;
end
 u=max(N);
mus=n1*u/n2;  
clear N n0 u;

C=rand(n2,1,n3);

X=c_product(Ur,C); % S子空间的信号
clear C;
C=rand(n1-n2,1,n3); %S的正交子空间的信号
Y=c_product(U_r,C);
clear C;

X0=0;
for j=1:n3
    a=norm(X(:,1,j));
    b=a^2;
    X0=X0+b;
end
Y0=0;
for j=1:n3
    a=norm(Y(:,1,j));
    b=a^2;
    Y0=Y0+b;
end
clear a b;

%A_c=c_mat(A);
P_c=c_mat(P);
X_u=unfold(X);
Y_u=unfold(Y);

E_E0=zeros(m,T);
Es_Es0=zeros(m,T);    

for t=1:T
    W=zeros(m,1);
    E=zeros(m,1);
    %E0=zeros(m,1);
    Es=zeros(m,1);
   
    for k=1:50:m
        W(1:k)=randperm(m,k);
        w=sort(W(1:k));
        Yw=zeros(k,1);
        Xw=zeros(k,1);
        
        Uw=zeros(k,n);
        
        for i=1:k
            Yw(i)=Y_u(w(i));
            Xw(i)=X_u(w(i));
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
        as=Pw*Xw;
        
        
        e=Yw-a;
        es=Xw-as;
                 
        clear a as  w Xw Yw Uw;
        %计算能量
        a=norm(e);
        E(k)=a^2;
       % E0(k)=norm(Yw);
       b=norm(es);
        Es(k)=b^2;
       % Es0(k)=norm(Xw);
      
        
      %  E1s0(k)=Es0(k);
      clear a b;
        E(k)
        %归一化
        E_E0(k,t)=E(k)/Y0;
        Es_Es0(k,t)=Es(k)/X0;
        
    end
end

Lmean=zeros(m,1);
Lmax=zeros(m,1);
Lmin=zeros(m,1);
Lsmean=zeros(m,1);
Lsmax=zeros(m,1);
Lsmin=zeros(m,1);

for i=1:m
    Lmean(i)=mean(E_E0(i,:));
    Lmax(i)=max(E_E0(i,:));
    Lmin(i)=min(E_E0(i,:));
    Lsmean(i)=mean(Es_Es0(i,:));
    Lsmax(i)=max(Es_Es0(i,:));
    Lsmin(i)=min(Es_Es0(i,:));
   
end

mu

k=1:50:m;
figure;
plot(k,Lmin(k),'x',k,Lmax(k),'+',k,Lmean(k),'*');axis([0 m 0 1.2]);xlabel('m'),ylabel('projection residual'),title('Projection Residual');
legend('minimum','maximum','mean');
figure;
plot(k,Lsmin(k),'x',k,Lsmax(k),'+',k,Lsmean(k),'*');axis([0 m 0 0.2]);xlabel('m'),ylabel('projection residual'),title('Projection Residual');
legend('minimum','maximum','mean');



