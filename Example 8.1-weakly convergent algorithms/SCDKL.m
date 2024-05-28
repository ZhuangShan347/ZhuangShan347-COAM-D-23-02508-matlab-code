function [iter,MSE,t,E,x1]=SCDKL(A,xs,b,x0,x1,kmax,k,pp)
E=(norm(x1-xs))/(max(1,norm(x1)));
iter=1;
MSE=(norm(x1-xs))^2/k;
while 1
    x=x1-x0;
    normx=norm(x);
    normx2=normx^2;
    m1=iter^2*max(normx2,normx);
    if x==0
        theta=0.27;
    else
        theta=min( 0.27,1/m1);
    end
    y=x1+theta*x;
    x0=x1;
    ab=A*y-b;
    grad_f=A'*ab; %% P_Q(y)=b,
    eta=max(1,norm(grad_f));
    lamda=0.1/eta;
    dd=lamda*grad_f;
    x1=y-dd; %%接下来计算在C上的投影
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
end
t=toc;
