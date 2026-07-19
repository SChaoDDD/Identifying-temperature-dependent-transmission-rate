function dy=seir(t,y,param)
global temp_ia
temp=temp_ia(floor(t+1));
beta0=param(5);
delta=0.0078/7;
sigma=1/2;
gamma=1/7;
dy=zeros(5,1);
S=y(1);
E=y(2);
I=y(3);
R=y(4);
N=12830000;
% k=0.0625;
% tempsize=1.3;
x1=temp;
beta=beta0*(exp((x1 .* (-0.155299307 .* exp(0.113332 .* x1))) ./ (exp(0.181309 .* x1) + 28.74496471)) ./ 0.804644);
beta=beta*(0.281824*cos(2*pi*t/365)+ 1.052856);
% if t>8
% if temp_ia(floor(t+7))-temp_ia(floor(t-7))<-tempsize 
%     beta=beta*(1+k);
% elseif temp_ia(floor(t+7))-temp_ia(floor(t-7))>tempsize
%     beta=beta*(1-k);
% end
% end
dy(1)=-beta*S*I/N+delta*R;
dy(2)=beta*S*I/N-sigma*E;
dy(3)=sigma*E-gamma*I;
dy(4)=gamma*I-delta*R;
dy(5)=sigma*E;
end