function [k,t,F,SUDI]=PC(A,x1,x0,lamda0)
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
    y=x1+theta*x;
    x0=x1;
    y1=A*y;   %%%AY
    yy=y1(2:20);
    Q=y1(1)+(yy'*yy)-1;%%%Q(AY)
    v=[1;2*yy]; %%% 计算出梯度
    %     normv=norm(v);
    v2=v'*v;
    if Q>0
        pq=y1-Q*v/v2;
    else              %%%%%%%%% compute AY project Qn；
        pq=y1;
    end
    %    f=(y1-pq)'*(y1-pq); %%去掉
    grad_fn=A'*(y1-pq);%% 计算函数在 y_n处的梯度 相当于Phi
    dd=lamda0* grad_fn;
    x=y-dd; %%% 接下来计算 C上的投影
    y2=y(2:10); %% n 在变化的时候，这里也在变化
    c=-y(1)+(y2'*y2);
    u=[-1; 2*y2];%%%grad_c
    u2=u'*u;
    c1=c-u'*dd; %% 出现了问题
    if c1>0
        x1=x-(c1*u)/u2;
    else
        x1=x;
    end
    y2=A*x1;
    yy=y2(2:20);
    Q=y2(1)+(yy'*yy)-1;%%%Q(AY)
    v=[1;2*yy]; %%% 计算出梯度
    %     normv=norm(v);
    v2=v'*v;
    if Q>0
        pq=y2-Q*v/v2;
    else              %%%%%%%%% compute AY project Qn；
        pq=y2;
    end
    grad_fh=A'*(y2-pq);%% 计算函数在 y_n处的梯度 相当于Phi
    f=0.499*norm(y-x1)^2;
    ff=dot(grad_fn-grad_fh,y-x1);
    p=(22/(k+1))^(1.1)+1;
    lamda=min(f/ff,lamda0*p);
    d=y-x1-lamda*(grad_fn-grad_fh);
    alp=dot(y-x1,d)+0.5*lamda*norm(y2-pq)^2;
    alp=alp/norm(d)^2;
    x1=y-0.5*alp*d;
    alpha=1/k;
    x1=alpha*x0+(1-alpha)*x1;
    if  norm(x1-y)<epss
        break;
    end
    k=k+1;
    F(k)=norm(x1-x0);
    SUDI(k)=norm(x1-y);
    lamda0=lamda;
end
t=toc;
end