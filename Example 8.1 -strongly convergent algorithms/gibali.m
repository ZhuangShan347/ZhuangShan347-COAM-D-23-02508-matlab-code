function [iter,MSE,t,SUDI,E,x1,F]=gibali(A,xs,b,x0,x1,kmax,k,pp)
E=(norm(x1-xs))/(max(1,norm(x1)));
iter=1;
MSE=(norm(x1-xs))^2/k;
SUDI(iter)=E;
while 1
    x=x1-x0;
    m1=norm(x)*(1/iter)^5;
    if x==0
        theta=0.9;
    else
        theta=min( 0.9,1/m1);
    end
    theta=max(0,theta-0.1);
    y=x1+theta*x;
    x0=x1;
    ab=A*y-b;
    f=(ab)'*(ab);
    grad_f=A'*ab; %% P_Q(y)=b,
    ksi=sign(y);
    ksi2=ksi'*ksi;
    c=norm(y,1)-pp;
    x1=y;
    if c>0
        x1=y-c*ksi/ksi2;
    end
    grad_g=y-x1;
    fg=norm(grad_f)^2+norm(grad_g)^2;
    lamda=2*f/fg;
%     eta=max(1,norm(grad_f));
%     lamda=0.1/eta;
    dd=lamda*grad_f;
    x1=y-dd; %%接下来计算在C上的投影
    ksi=sign(y);
    ksi2=ksi'*ksi;
    c=norm(y,1)-pp-ksi'*dd;
    if c>0
        x1=x1-c*ksi/ksi2;
    end
    alpha=0.00001/iter;
    beta=0.7;
    x1=(1-alpha-beta)*y+beta*x1;
    E=(norm(x1-xs))/(max(1,norm(x1)));
    MSE=(norm(x1-xs)^2)/k;
    if iter>kmax
        break;
    end
    iter=iter+1;
    SUDI(iter)=E;
    F(iter)=lamda;
end
t=toc;
