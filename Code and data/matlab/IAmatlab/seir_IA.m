function dy=seir_IA(t,y,param)
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
N=3125800;
x1=temp;
k=1.1005;
temperature=temp;
%beta=beta0*exp((((x1 .* -0.0006722514891077485) + -0.008756861382359383) .* x1) + 0.22189917446131707);%up

beta=beta0*(exp((x1 .* (-0.155299307 .* exp(0.113332 .* x1))) ./ (exp(0.181309 .* x1) + 28.74496471)) ./ 0.804644);%mean
beta=beta*(0.291786*cos(2*pi*t/365)+1.006805);
% if t>10
% if mean(temp_ia(floor(t+3):floor(t+6)))-mean(temp_ia(floor(t-6):floor(t-3)))<-tempsize
%     beta=beta*k;
% elseif mean(temp_ia(floor(t+3):floor(t+6)))-mean(temp_ia(floor(t-6):floor(t-3)))>tempsize
%     beta=beta*(2-k);
% end
% end
dy(1)=-beta*S*I/N+delta*R;
dy(2)=beta*S*I/N-sigma*E;
dy(3)=sigma*E-gamma*I;
dy(4)=gamma*I-delta*R;
dy(5)=sigma*E;
end