function [ c ] = c_mat( a)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明
[n1 n2 n3]=size(a);
m=n1*n3; n=n2*n3;
A1=zeros(m,n); % 存储 Toeplitz 矩阵
A2=zeros(m,n);  % 存储 Hankel 矩阵

for j=1:n3           %计算 Toeplitz 矩阵
    for i=1:n3
        if i==j
            A1((i-1)*n1+1:i*n1,(j-1)*n2+1:j*n2)=a(:,:,1);
        end
        if i<j
           A1((i-1)*n1+1:i*n1,(j-1)*n2+1:j*n2)=a(:,:,j-i+1);
        end
        if i>j
          A1((i-1)*n1+1:i*n1,(j-1)*n2+1:j*n2)=a(:,:,i-j+1);
        end
        
    end  
end

for j=1:n3         %计算 Hankel 矩阵
    for i=1:n3
       k=i+j;
       if k<n3+1
           A2((i-1)*n1+1:i*n1,(j-1)*n2+1:j*n2)=a(:,:,k);
       end
       if k==n3+1
           A2((i-1)*n1+1:i*n1,(j-1)*n2+1:j*n2)=zeros(n1,n2);
       end
       if k>n3+1
           k=2*(n3+1)-k;
          A2((i-1)*n1+1:i*n1,(j-1)*n2+1:j*n2)=a(:,:,k);    
       end
        
    end  
end

b=A1+A2;  %输出 Toeplitz 矩阵 + Hankel 矩阵

I_n3=eye(n3);
I_n1=eye(n1);
I_n2=eye(n2);
Z=diag(ones(n3-1, 1), 1);
IZ=I_n3+Z;
IZ_INV=inv(IZ);
F=kron(IZ_INV,I_n1);
B=kron(IZ,I_n2);

c=F*b*B;

end

