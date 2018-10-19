function [ b ] = c_transpos( a )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
[n1 n2 n3]=size(a);
b=zeros(n1,n2,n3);
A=zeros(n1,n2,n3);
B=zeros(n2,n1,n3);
A=my_dct(a);

     for k=1:n3
        B(:,:,k)=A(:,:,k)'; 
     end
     
     b=my_idct(B);

end

