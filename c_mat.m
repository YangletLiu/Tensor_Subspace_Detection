function [ c ] = c_mat( a)
%UNTITLED8 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[n1 n2 n3]=size(a);
m=n1*n3; n=n2*n3;
A1=zeros(m,n); % �洢 Toeplitz ����
A2=zeros(m,n);  % �洢 Hankel ����

for j=1:n3           %���� Toeplitz ����
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

for j=1:n3         %���� Hankel ����
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

b=A1+A2;  %��� Toeplitz ���� + Hankel ����

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

