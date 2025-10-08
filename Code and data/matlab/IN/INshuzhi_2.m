clear 
clc
global temp_ia_1 I_IN
data1=readtable('IN.csv');
data1=table2array(data1(1:308,4));
I_IN=data1;
data2=readtable('INclimate.xlsx');
data2=table2array(data2(273:273+2200,3));
temp_in=data2(1:2200);
window_size = 7;

temp_ia_1 = zeros(length(temp_in) , 1);
for i = window_size:length(temp_ia_1)
    temp_ia_1(i) = mean(temp_in(i-window_size+1:i));
end
temp_ia_1=temp_ia_1(window_size:end);
N=6570000;
E0=10;
I0=10;
R0=9;
S0=N-E0-I0-R0;
beta0=0.1455;
delta=0.0078/7;
sigma=1/2;
gamma=1/7;
y0=[S0,E0,I0,R0,0];
x0=[S0,E0,I0,R0,beta0];
%x0=[S0,E0,I0,R0,beta0,k0];
% options=optimset('Display','on','MaxIter', 5000, 'MaxFunEvals',3000);
% [x,fval]=fmincon(@fun_IN,x0,[],[],[],[],low,up,[],options);
low=[1,1,1,1,0];
up=[3094000,2000,3500,5000,1];
tend=307*7;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
sol=ode45(@seir_IN,[0 tend-1],y0,options,x0);
x1=linspace(0,tend-1,307);
y1=deval(sol,x1);
newcase=y1(5,3:end)-y1(5,2:end-1);
procase=load('notempcase.mat');
procase=procase.newcase;
week1 = 1:length(newcase);
week2 = 1:307;
startYear = 2013;
startWeek = 40;
startDate = datetime(startYear, 1, 1) + calweeks(startWeek - 1);
date1 = startDate + calweeks(week1 - 1);
date2 = startDate + calweeks(week2 - 1);
figure;
plot(date1,procase,'LineWidth',1.5,'LineStyle','--','Color','k')
hold on
plot(date1,newcase,'LineWidth',1.5,'Color','b')
hold on
plot(date2,I_IN(1:307),'r.', 'MarkerSize', 15)
xtickformat('yyyy')
% xtickangle(45);
xlabel('year');
ylabel('case');
legend('unadjust','adjust','case')




