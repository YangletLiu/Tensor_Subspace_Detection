function [ a] = dct_trans( y)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
         n = max(size(y));
         x = y(:);
         a=zeros(1,1,n);

         I = eye(n);
         C = dct(eye(n));
         Z = diag(ones(n-1, 1), 1);
         W = diag(C(:, 1));

         b = inv(W) * C * (I + Z) * x;
         
         for i=1:n
             a(1,1,i)=b(i);
         end
end

