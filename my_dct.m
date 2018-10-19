function [ A ] = my_dct( a )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    [n1 n2 n3]=size(a);
     A=zeros(n1,n2,n3);
     for i=1:n1
           for j=1:n2
               g=dct_trans(a(i,j,:));
               A(i,j,:)=g;
           end
      end
       clear g;

end

