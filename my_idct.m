function [ a ] = my_idct( A )
%UNTITLED4 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    [n1 n2 n3]=size(A);
     a=zeros(n1,n2,n3);
     for i=1:n1
           for j=1:n2
               g=idct_trans(A(i,j,:));
               a(i,j,:)=g;
           end
      end
       clear g;

end

