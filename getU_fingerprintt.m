clear;

T_data=importdata('tensor fingerprint data\tensor.mat');
V_data=unfold(T_data);

n1=64;n2=256;n3=21;
threshold=1e-8;
train_num=10;

Energy=zeros(n2,1);
for i=1:n2
    Energy(i)=norm(V_data(:,i));
    T_data(:,i,:)=T_data(:,i,:)/Energy(i);
    V_data(:,i,:)=V_data(:,i,:)/Energy(i);
end


V_size=n1*n3;

t_train=T_data(:,1:train_num,:);

%处理训练数据
[T_U, T_S, T_V]=t_svd(t_train);


S_E=zeros(train_num,1);

for i=1:train_num
    if i<=n1
        E=0;
        for j=1:n3
        E=E+T_S(i,i,j)^2;
        end
        S_E(i)=E^2;
    end
end

t_r=9;

t_U=T_U(:,1:t_r,:);
t_S=T_S(1:t_r,1:t_r,:);
t_V=T_V(:,1:t_r,:);
a=transpos(t_V);
b=t_product(t_U,t_S);
T_train=t_product(b,a);

 V_train=unfold(T_train);
 
 [V_U,V_S,V_V]=svd(V_train);
 v_r=10;
 v_U=V_U(:,1:v_r);
 
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 index=1:n2;
index_s_=zeros(1,n2);

s_num=0;

 a=transpos(t_U);
 b=t_product(a,t_U);
 f=fft(b,[],3);
 g=zeros(t_r,t_r,n3);
 for i=1:n3
      g(:,:,i)=pinv(f(:,:,i));
 end
 h=ifft(g,[],3);
 s=t_product(t_U,h);
 t_Ps=t_product(s,a);
 clear a b f g h s i;
 
 %计算v_Ps
 a=v_U';
 b=a*v_U;
 f=pinv(b);
 g=v_U*f;
 v_Ps=g*a;
 clear a b f g v;

 T_E=zeros(n2,1);

 
for i=1:n2
    t_e=T_data(:,i,:)-t_product(t_Ps,T_data(:,i,:));
    unfold_te=unfold(t_e);
    a=norm(unfold_te);
    T_E(i)=a^2;
    if T_E(i)<threshold;
        s_num=s_num+1;
        index_s_(s_num)=i;
        
    end
end

index_s=index_s_(1,1:s_num);
index_ors=setdiff(index,index_s);

T_data_s=zeros(n1,s_num,n3);
T_data_ors=zeros(n1,n2-s_num,n3);
for i=1:s_num
    T_data_s(:,i,:)=T_data(:,index_s(i),:);
end
ors_num=n2-s_num;
for i=1:ors_num
   T_data_ors(:,i,:)=T_data(:,index_ors(i),:);
end
V_data_s=unfold(T_data_s);
V_data_ors=unfold(T_data_ors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 V_E=zeros(n2,1);
 v_e=V_data-v_Ps*V_data;
  
  for i=1:n2
      E=norm(v_e(:,i));
      V_E(i)=E^2;
   end


