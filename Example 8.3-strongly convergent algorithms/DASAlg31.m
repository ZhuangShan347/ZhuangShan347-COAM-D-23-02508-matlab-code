function [k,t,F,SUDI]=DASAlg31(A,x1,x0)
clc
global epss
tic
k=1;
theta=0.3;%% 0.5
beta=5;%% 10
lamda=20; %% 14
t=1;
xs=x0;
while 1
    y=x1+theta*(x1-x0);
    x0=x1;
    y2=y(2:10); %% n 在变化的时候，这里也在变化
    c=-y(1)+(y2'*y2);
    u=[-1; 2*y2];%%%grad_c
    u2=u'*u;
    c1=c; %% 出现了问题
    if c1>0
        w=y-(c1*u)/u2;
    else
        w=y;
    end
    y1=A*w;   %%%AY
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
    grad_fn=A'*(y1-pq);%% 计算函数在 y_n处的梯度
    dd=beta* grad_fn;
    x=w-dd; %%% 接下来计算 C上的投影
    y2=y(2:10); %% n 在变化的时候，这里也在变化
    c=-y(1)+(y2'*y2);
    u=[-1; 2*y2];%%%grad_c
    u2=u'*u;
    c1=c+u'*(x-y); %% 出现了问题
    if c1>0
        vv=x-(c1*u)/u2;
    else
        vv=x;
    end
    y1=A*vv;   %%%AY
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
    grad_v=A'*(y1-pq);%% 计算函数在 y_n处的梯度
    dist=dot( grad_fn,w-vv)-lamda*dot(grad_fn- grad_v,w-vv);
    while dist<0
        beta=0.01*beta;
        dd=beta* grad_fn;
        x=w-dd; %%% 接下来计算 C上的投影
        y2=y(2:10); %% n 在变化的时候，这里也在变化
        c=-y(1)+(y2'*y2);
        u=[-1; 2*y2];%%%grad_c
        u2=u'*u;
        c1=c+u'*(x-y); %% 出现了问题
        if c1>0
            vv=x-(c1*u)/u2;
        else
            vv=x;
        end
        y1=A*vv;   %%%AY
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
        grad_v=A'*(y1-pq);%% 计算函数在 y_n处的梯度
        dist=dot(grad_fn,w-vv)-lamda*dot(grad_fn- grad_v,w-vv);
    end
    z=w-t*dot( grad_v,w-vv)*grad_v/norm(grad_v)^2;
    y2=y(2:10); %% n 在变化的时候，这里也在变化
    c=-y(1)+(y2'*y2);
    u=[-1; 2*y2];%%%grad_c
    u2=u'*u;
    c1=c+u'*(z-y); %% 出现了问题
    if c1>0
        u=z-(c1*u)/u2;
    else
        u=z;
    end
    Ph=xs-max(0,xs-x1);
    d=0.5*norm(w-u)^2-dot(u-w,u);
    if dot(w-u,Ph)>d
        Ps=Ph-1.9*(dot(w-u,Ph)-d)*(w-u)/norm(w-u)^2;
    else
        Ps=Ph;
    end
    y2=y(2:10); %% n 在变化的时候，这里也在变化
    c=-y(1)+(y2'*y2);
    u=[-1; 2*y2];%%%grad_c
    u2=u'*u;
    c1=c+u'*(Ps-y); %% 出现了问题
    if c1>0
        x1=Ps-(c1*u)/u2;
    else
        x1=Ps;
    end
    if  norm(x1-y)<epss
        break;
    end
    k=k+1;
    F(k)=norm(x1-x0);
    SUDI(k)=norm(x1-y);
end
t=toc;
end