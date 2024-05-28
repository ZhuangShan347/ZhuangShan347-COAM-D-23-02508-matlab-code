function [k,t,F,SUDI]=gibali(A,x1,x0)
clc
global epss
tic
k=1;
while 1
    x=x1-x0;
    mm=norm(x)*(1/k)^3;
    if x~=0 
    theta=min(0.3,1/mm);
    end
    theta=max(0,theta-0.1);
    y=x1+theta*x;
    x0=x1;
    y1=A*y;   %%%AY
    yy=y1(2:20);
    Q=y1(1)+(yy'*yy)-1;%%%Q(AY)
    v=[1;2*yy]; %%% ������ݶ�
%     normv=norm(v);
    v2=v'*v;
    if Q>0
        pq=y1-Q*v/v2;
    else              %%%%%%%%% compute AY project Qn��
        pq=y1;
    end
    f=(y1-pq)'*(y1-pq);
    grad_fn=A'*(y1-pq);%% ���㺯���� y_n�����ݶ�
    ff=norm(grad_fn)^2;
    y2=y(2:10); %% n �ڱ仯��ʱ������Ҳ�ڱ仯
    c=-y(1)+(y2'*y2);
    u=[-1; 2*y2];%%%grad_c
    u2=u'*u;
    c1=c; %% ����������
    if c1>0
        x1=y-(c1*u)/u2;
    else
        x1=y;
    end
    grad_g=y-x1;
    fg=ff+norm(grad_g)^2;
    lamda=0.09*f/fg;
    dd=lamda* grad_fn;
    x=y-dd; %%% ���������� C�ϵ�ͶӰ
    y2=y(2:10); %% n �ڱ仯��ʱ������Ҳ�ڱ仯
    c=-y(1)+(y2'*y2);
    u=[-1; 2*y2];%%%grad_c
    u2=u'*u;
    c1=c-u'*dd; %% ����������
    if c1>0
        x1=x-(c1*u)/u2;
    else
        x1=x;
    end
     alpha=0.0001/k;
     beta=0.8;
     x1=(1-beta-alpha)*y+beta*x1;
     if  norm(x1-y)<epss
        break;
     end
      k=k+1; 
      F(k)=norm(x1-x0);
      SUDI(k)=norm(x1-y);
end 
t=toc;
end