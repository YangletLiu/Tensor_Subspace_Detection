tic;

n1=50; n2=10; n3=50; T=1; %参数
m=n1*n3;n=n2*n3;
%
A=rand(n1,n2,n3);
[U S V]=t_svd(A);
% X in S
Ur=U(:,1:n2,:);
c=rand(n2,1,n3);
X=t_product(Ur,c);
clear c;
X0=0;
for j=1:n3
    a=norm(X(:,1,j));
    b=a^2;
    X0=X0+b;
end
%正交投影
a=transpos(Ur);
Ps=t_product(Ur,a);
clear a;
%计算mu
N=zeros(n1,1);
for j=1:n1
n_=0;
for i=1:n1
a=0;
for k=1:n3
h=Ps(i,j,k)^2;
a=a+h;
end
n_=n_+a;
end
N(j)=n_;
end
c=max(N);
mu=n1*c/n2;
clear c n_ h a N;

% 正交子空间中的元素
U_r=U(:,n2+1:n1,:);
c=rand(n1-n2,1,n3);
Y=t_product(U_r,c);
clear c;
Y0=0;
for j=1:n3
    a=norm(Y(:,1,j));
    b=a^2;
    Y0=Y0+b;
end
clear a b;

A_c=bcirc(A);
P_c=bcirc(Ps);
X_u=unfold(X);
Y_u=unfold(Y);

E_E0=zeros(m,T);
Es_Es0=zeros(m,T);    
for t=1:T
    W=zeros(m,1);
    E=zeros(m,1);
    %E0=zeros(m,1);
    Es=zeros(m,1);
   % Es0=zeros(m,1);
    
    for k=100:100
       
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
      clear a b c;
        E(k)
        %归一化
        E_E0(k,t)=E(k)/Y0;
        Es_Es0(k,t)=Es(k)/X0;
       
        toc;
        
        disp(['运行时间: ',num2str(toc)]);
    end
end
