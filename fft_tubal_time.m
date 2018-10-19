n1=50; n2=10; n3=50; T=1;
A=rand(n1,n2,n3);
[U S V]=t_svd(A);
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
clear a b;
a=transpos(Ur);
Ps=t_product(Ur,a);
clear a;
U_r=U(:,n2+1:n1,:);
c=rand(n1-n2,1,n3);
y=t_product(U_r,c);
clear c;

Y0=0;
for j=1:n3
    a=norm(y(:,1,j));
    b=a^2;
    Y0=Y0+b;
end
clear a b;

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
Yw(i,1,:)=y(w(i),1,:);
Xw(i,1,:)=X(w(i),1,:);
Uw(i,:,:)=A(w(i),:,:);

end
tic;
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
toc;
tic;
clear a b f g h s;
a=t_product(Pw,Yw);
as=t_product(Pw,Xw);           
e=Yw-a;
es=Xw-as;

clear a as a1 a1s;
n=0;

ns=0;


for i=1:n3
a=norm(e(:,1,i));
as=norm(es(:,1,i));

b=a^2;
bs=as^2;

n=n+b;

ns=ns+bs;

end
toc;
E(k)=n;

Es(k)=ns;

E(k);
clear a a0 b b0 n n0 e as as0 bs bs0 ns ns0 es a1 a10 b1 b10 n_1 n10 e1 a1s a1s0 b1s b1s0 n1s n1s0 e1s;
 X1w=zeros(n1,1,n3);     %用来存储X采样后的数据（zero-filled)
 Y1w=zeros(n1,1,n3);     %用来存储Y采样后的数据（zero-filled)
end

E_Y0(:,t)=E./Y0;
Es_X0(:,t)=Es./X0;    
        %
end
