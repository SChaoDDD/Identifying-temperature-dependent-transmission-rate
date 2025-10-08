function dy=seir_SD(t,y,param)
global temp_ia_1
temp=temp_ia_1(floor(t+1));
beta0=param(5);
delta=0.0078/7;
sigma=1/2;
gamma=1/7;
b=0;
%c=param(10);
dy=zeros(5,1);
S=y(1);
E=y(2);
I=y(3);
R=y(4);
N=866000;
x=temp;
k=1.1236;
beta=beta0*(exp((-8.96*10^(-7).*x.^3-2.35*10^(-4).*x.^2-6.03*10^(-3).*x)+0.1223196670929));
if temp_ia_1(floor(t+30))-temp_ia_1(floor(t+1))<-3
    beta=beta*k;
elseif temp_ia_1(floor(t+30))-temp_ia_1(floor(t+1))>3
    beta=beta*(2-k);
end
%beta=beta0*(ppval(pp, temp));
dy(1)=-beta*S*I/N+delta*R;
dy(2)=beta*S*I/N-sigma*E;
dy(3)=sigma*E-gamma*I;
dy(4)=gamma*I-delta*R;
dy(5)=sigma*E;
end