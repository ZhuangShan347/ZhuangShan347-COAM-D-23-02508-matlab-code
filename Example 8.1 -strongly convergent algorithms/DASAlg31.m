function [iter,MSE,t,SUDI,E,x1]=DASAlg31(A,xs,b,x0,x1,kmax,k,pp)
E=(norm(x1-xs))/(max(1,norm(x1)));
iter=1;
MSE=(norm(x1-xs))^2/k;
SUDI(iter)=E;
theta=0.7;%% 0.5
beta=8;%% 10
lamda=12; %% 14
while 1
    x=x1-x0;
    y=x1+theta*x;
    x0=x1;
    ksi=sign(y);
    ksi2=ksi'*ksi;
    c=norm(y,1)-pp;%%接下来计算y处的投影
    if c>0
        w=y-c*ksi/ksi2;
    else 
        w=y;
    end
    ab=A*w-b;
    grad_f=A'*ab; %% P_Q(y)=b,grad_f=Fw
    z=w-beta* grad_f;
    c=norm(y,1)-pp+ksi'*(z-y);%%接下来计算z处的投影
    if c>0
        v=z-c*ksi/ksi2;
    else
        v=z;
    end
   Fv=A'*(A*v-b);
    xx=dot(grad_f,w-v)-lamda*dot(grad_f-Fv,w-v);
    while xx<0
        beta=0.9*beta; %% 0.6 0.9
        z=w-beta* grad_f;
    c=norm(y,1)-pp+ksi'*(z-y);%%接下来计算z处的投影
    if c>0
        v=z-c*ksi/ksi2;
    else
        v=z;
    end
   Fv=A'*(A*v-b);
    xx=dot(grad_f,w-v)-lamda*dot(grad_f-Fv,w-v);
    end 
      zz=w-dot(Fv,w-v)*Fv/norm(Fv)^2;
      c=norm(y,1)-pp+ksi'*(zz-y);%%接下来计算z处的投影
    if c>0
        u=zz-c*ksi/ksi2;
    else
        u=zz;
    end
    Ph=xs-max(0,xs-x1); %% 在 H上投影   
    d=0.5*norm(w-u)^2-dot(u-w,u);
    if dot(w-u,Ph)>d
        Ps=Ph-1.9*(dot(w-u,Ph)-d)*(w-u)/norm(w-u)^2;
    else
        Ps=Ph;
    end
     c=norm(y,1)-pp+ksi'*(Ps-y);
   
    if c>0
        x1=Ps-c*ksi/ksi2;
    else
        x1=Ps;
    end
    E=(norm(x1-xs))/(max(1,norm(x1)));
    MSE=(norm(x1-xs)^2)/k;
    if iter>kmax
        break;
    end
    iter=iter+1;
    SUDI(iter)=E;
end
t=toc;
end