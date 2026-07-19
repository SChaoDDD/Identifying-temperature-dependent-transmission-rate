clear 
clc
global temp_ia
data1=readtable('IL.csv');
timemark=574;
data1=table2array(data1(1:timemark-1,6));
I_il=data1.*2.74;%每周流感数据
data3=readtable('ILaverage8.csv');
data3=table2array(data3(:,2));
temp_ia=data3;
% data2=readtable('ILdaily_beta_ensemble1.csv');
% beta1=table2array(data2(:,3));

% yyaxis left
% plot(ppval(pp,data3(1:2184)))
% yyaxis right
% plot(data3(1:2184))

N=12830000;

E0=8;
I0=6;
R0=0;
S0=N-E0-I0-R0;
beta0=0.1434;%0.1485
delta=0.0078/7;
sigma=1/2;
gamma=1/7;
y0=[S0,E0,I0,R0,0];
x0=[S0,E0,I0,R0,beta0];
low=[1,10,10,50,0];
up=[12900000,2000,3500,5000,1];%先给一个上下界然后后面再调
tend=timemark*7;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
sol=ode45(@seir,[0 tend-1],y0,options,x0);
%valid_tspan = tspan(tspan >= min(sol.x) & tspan <= max(sol.x));
x1=linspace(0,tend-1,timemark); %按周为一个时间单位
y1=deval(sol,x1);
% plot(y1(3,:))
% hold on
% plot(y1(2,:))
newcase=y1(5,2:end)-y1(5,1:end-1);
week1 = 1:timemark-1;
week2 = 1:timemark-1;
startYear = 2013;
startWeek = 40;
startDate = datetime(startYear, 1, 1) + calweeks(startWeek - 1);
date1 = startDate + calweeks(week1 - 1);
date2 = startDate + calweeks(week2 - 1);
figure;
plot(date1,newcase,'LineWidth',1.5,'Color','b')
hold on
plot(date2,I_il(1:timemark-1),'r.', 'MarkerSize', 15)
xtickformat('yyyy')
xlabel('year');
ylabel('ILI cases');
legend('Unadjusted curve','Adjusted curve','Real data')



% IAcase=[I_il(1:307),newcase'];
save('ILfit.mat',"newcase")

% indown=load("ILdowncase.mat");
% inmean=load("ILmeancase.mat");
% inup=load("ILroiginalcase.mat");
% indown=indown.IAcase;
% inmean=inmean.IAcase;
% inup=inup.IAcase;
% figure
% hold on
% plot(date2,inmean(1:307,1),'.','Color','#E63946','MarkerSize',18);       % 大红散点
% plot(date2,inmean(1:307,2),'-','Color','#0066FF','LineWidth',2);  % 深蓝均值线
% plot(date2,inup(1:307,2),'--','Color','#0066FF','LineWidth',2);
% % x = 1:length(inmean(1:307,2));
% % y_low  = indown(1:307,2);
% % y_high = inup(1:307,2);
% % 
% % %亮黄色阴影
% % fill([date2, fliplr(date2)], [y_low; flipud(y_high)], ...
% %     [0.2, 0.6, 1], 'FaceAlpha', 0.45, 'EdgeColor','none');
% 
% hold off
% grid on
% xtickformat('yyyy')


