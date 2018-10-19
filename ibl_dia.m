function [ a ] = ibl_dia( A,n1,n2,n3 )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
     a=zeros(n1,n2,n3);
     for k=1:n3
         a(:,:,k)=A((k-1)*n1+1:k*n1,(k-1)*n2+1:k*n2);
     end


end

