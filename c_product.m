function [ c ] = c_product( a,b )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[a_n1 a_n2 a_n3]=size(a);
    [b_n1 b_n2 b_n3]=size(b);
    if a_n2==b_n1&&a_n3==b_n3
       
       A=my_dct(a);
       B=my_dct(b);
       C=zeros(a_n1,b_n2,a_n3);
       
      for k=1:a_n3
          C(:,:,k)=A(:,:,k)*B(:,:,k);
       end
       
      c=my_idct(C);
      
    else
        c=0;
    end

end

