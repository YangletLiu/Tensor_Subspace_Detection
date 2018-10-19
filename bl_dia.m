function [ A ] = bl_dia( a )
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
     [n1 n2 n3]=size(a);
     m=n1*n3; n=n2*n3;
     A=zeros(m,n);
     for k=1:n3
        A((k-1)*n1+1:k*n1,(k-1)*n2+1:k*n2)=a(:,:,k); 
     end

end

