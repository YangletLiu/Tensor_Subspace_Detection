
clear;

load data_fingerprint;

S_N=1:64;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%¼ì²â¸ÅÂÊ

t_Re_ors=zeros(1,64);
v_Re_ors=zeros(1,64);

for ii=1:64
    K=S_N(ii);
    Th=threshold*K/n1;
        
    num=0;
    counter1=0;
    counter2=0;
    while (counter1<100&&counter2<100)&&num<1000

          W(1:K)=randperm(n1,K);
          w=sort(W(1:K));
          T_ors_w=zeros(K,ors_num,n3);
          t_Uw=zeros(K,t_r,n3);
          V_ors_w=zeros(K*n3,ors_num);
          v_Uw=zeros(K*n3,v_r);
          for i=1:K
              T_ors_w(i,:,:)=T_data_ors(w(i),:,:);
              t_Uw(i,:,:)=t_U(w(i),1:t_r,:);
              for j=1:n3
                  V_ors_w(i+(j-1)*K,:)=V_data_ors(w(i)+(j-1)*K,:);
                  v_Uw(i+(j-1)*K,:)=v_U(w(i)+(j-1)*K,1:v_r);
              end
           end
           clear i j;
           %¼ÆËãt_Pw
           a=transpos(t_Uw);
           b=t_product(a,t_Uw);
           f=fft(b,[],3);
           g=zeros(t_r,t_r,n3);
           for i=1:n3
               g(:,:,i)=pinv(f(:,:,i));
           end
           h=ifft(g,[],3);
           s=t_product(t_Uw,h);
           t_Pw=t_product(s,a);
           clear a b f g h s i;

            %¼ÆËãv_Pw
            a=v_Uw';
            b=a*v_Uw;
            f=pinv(b);
            g=v_Uw*f;
            v_Pw=g*a;
            clear a b f g v;

            t_e=T_ors_w-t_product(t_Pw,T_ors_w);
            v_e=V_ors_w-v_Pw*V_ors_w;

            unfold_te=unfold(t_e);

            T_E=zeros(s_num,1);
            V_E=zeros(s_num,1);
            for i=1:s_num
                E=norm(unfold_te(:,i));
                T_E(i)=E^2;
                if T_E(i)>=Th
                    counter1=counter1+1;
                end
                E=norm(v_e(:,i));
                V_E(i)=E^2;
                if V_E(i)>=Th
                    counter2=counter2+1;
                end
                num=num+1;
            end
            counter1
            counter2
            num

        end;
        t_Pd=counter1/num;
        v_Pd=counter2/num;
        t_Re_ors(ii)=t_Pd;
        v_Re_ors(ii)=v_Pd;

    
    
end

i=1:4:64;
figure;

plot(i,t_Re_ors(i),'-o');xlabel('Sampling Number');ylabel('P_{D}');title('Detection Probability');
hold on;
plot(i,v_Re_ors(i),'-v');
axis([0 64 0 1]);
legend('tensor','vector');
grid minor;