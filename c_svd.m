function [ u s v ] = c_svd( a)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    [n1 n2 n3]=size(a);
    u=zeros(n1,n1,n3);
    s=zeros(n1,n2,n3);
    v=zeros(n2,n2,n3);
    U=zeros(n1,n1,n3);
    S=zeros(n1,n2,n3);
    V=zeros(n2,n2,n3);
    A=zeros(n1,n2,n3);
    
    A=my_dct(a);
    
      for k=1:n3
       [b c d]=svd(A(:,:,k));
       U(:,:,k)=b;
       S(:,:,k)=c;
       V(:,:,k)=d;
       clear b c d;
      end
       
      u=my_idct(U);
      s=my_idct(S);
      v=my_idct(V);

end

