function [iter,MSE,t,SUDI,E,x1]=SYGAAlg1(A,xs,b,x0,x1,kmax,k,pp,lamda0)
E=(norm(x1-xs))/(max(1,norm(x1)));
iter=1;
MSE=(norm(x1-xs))^2/k;
SUDI(iter)=E;
delta=0.5;
while 1
    x=x1-x0;
    if mod(iter,2)==0
        y=x1;
    else
      y=x1+0.3*x;
    end
    x0=x1;
    ab=A*y-b;
    grad_f=A'*ab; %% P_Q(y)=b,
    dd=lamda0*grad_f;
    z=y-dd; %%接下来计算在C上的投影
    ksi=sign(y);
    ksi2=ksi'*ksi;
    c=norm(y,1)-pp-ksi'*dd;
    if c>0
       z=z-c*ksi/ksi2;
    end
    xx=lamda0*norm(A'*(A*y-A*z))-delta*norm(y-z);
    while xx>0
        lamda0=0.5*lamda0;
         dd=lamda0*grad_f;
    z=y-dd; %%接下来计算在C上的投影
    ksi=sign(y);
    ksi2=ksi'*ksi;
    c=norm(y,1)-pp-ksi'*dd;
    if c>0
       z=z-c*ksi/ksi2;
    end
    xx=lamda0*norm(A'*(A*y-A*z))-delta*norm(y-z);
    end
    ab=A*z-b;
    grad_f=A'*ab; %% P_Q(y)=b,
    dd=lamda0*grad_f;
    x1=y-dd;
    ksi=sign(y);
    ksi2=ksi'*ksi;
    c=norm(y,1)-pp-ksi'*dd;
    if c>0
       x1=x1-c*ksi/ksi2;
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
