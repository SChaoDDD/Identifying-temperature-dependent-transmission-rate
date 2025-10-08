clear 
clc
global temp_ia_1 I_il
data1=readtable('IAliugan.csv');
data1=table2array(data1);
data1=data1(1:end,:);
I_il=data1(:,5);
data2=readtable('climate_IA.csv');
data2=table2array(data2);
temp_ia=data2(1:end,2);
window_size = 7;

temp_ia_1 = zeros(length(temp_ia) , 1);
for i = window_size:length(temp_ia_1)
    temp_ia_1(i) = mean(temp_ia(i-window_size+1:i));
end
temp_ia_1=temp_ia_1(window_size:end);
N=3094000;
E0=5;
I0=10;
R0=10;
S0=N-E0-I0-R0;
beta0=0.1450;%0.1509
delta=0.0078/7;
sigma=1/2;
gamma=1/7;
b=0;

y0=[S0,E0,I0,R0,0];
x0=[S0,E0,I0,R0,beta0];
low=[1,1,1,1,0];
up=[3094000,2000,3500,5000,1];
tend=582*7;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
sol=ode45(@seir_IA,[0 tend-1],y0,options,x0);
x1=linspace(0,tend-1,582); 
y1=deval(sol,x1);
newcase=y1(5,3:end)-y1(5,2:end-1);
week1 = 1:305;
week2=305:580;
week3 = 1:580;
startYear = 2013;
startWeek = 40;
startDate = datetime(startYear, 1, 1) + calweeks(startWeek - 1);
date1 = startDate + calweeks(week1 - 1);
date2 = startDate + calweeks(week2 - 1);
date3 = startDate + calweeks(week3 - 1);
figure;
plot(date1,newcase(1:305),'LineWidth',1.5,'Color','k')
hold on
plot(date2,newcase(305:580),'LineWidth',1.5,'Color','b')
hold on
plot(date3,I_il(1:580),'r.', 'MarkerSize', 15)
xtickformat('yyyy')
% xtickangle(45);
xlabel("year");
ylabel('case');






