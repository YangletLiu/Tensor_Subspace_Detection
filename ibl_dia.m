function [ a ] = ibl_dia( A,n1,n2,n3 )
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
     a=zeros(n1,n2,n3);
     for k=1:n3
         a(:,:,k)=A((k-1)*n1+1:k*n1,(k-1)*n2+1:k*n2);
     end


end

