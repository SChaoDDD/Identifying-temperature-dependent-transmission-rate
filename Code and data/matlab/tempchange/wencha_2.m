clear 
clc

% iatemp=readtable("climate_IA.csv");
% iatemp=table2array(iatemp(:,2));
% weektemp=zeros(308,1);
% for i=1:308
%     for j=(i-1)*7+1:i*7
%         weektemp(i)=weektemp(i)+iatemp(j);
%     end
% end
% weektemp=weektemp./7;
% d_iatemp=weektemp(2:end)-weektemp(1:end-1);
% beta=readtable("shuju.csv");
% iabeta=table2array(beta(:,6));
% iabeta=iabeta(2:end);
% temp_beta_1=zeros(307,2);
% temp_beta_2=zeros(307,2);
% averbeta=readtable("IAbeta.csv");

iatemp=readtable("INclimate.xlsx");
iatemp=table2array(iatemp(273:273+2200,3));
weektemp=zeros(308,1);
for i=1:308
    for j=(i-1)*7+1:i*7
        weektemp(i)=weektemp(i)+iatemp(j);
    end
end
weektemp=weektemp./7;
d_iatemp=weektemp(2:end)-weektemp(1:end-1);

beta=readtable("INbeta.xlsx");
iabeta=table2array(beta(:,1));
temp_beta_1=zeros(307,2);
temp_beta_2=zeros(307,2);
iabeta=iabeta(2:end);
temp_beta=zeros(307,2);
averbeta=readtable("INquxian.csv");
x0=table2array(averbeta(:,1));
y0=table2array(averbeta(:,2));
pp = spline(x0, y0);
% temp_beta(:,1)=weektemp(1:307);
% temp_beta(:,2)=iabeta(1:307);
for i=1:303
    if weektemp(i)-weektemp(i+4)>-2.5&&weektemp(i)-weektemp(i+4)<2.5
        temp_beta_1(i,1)=weektemp(i);
        temp_beta_1(i,2)=iabeta(i)/0.1455*7*ppval(pp,weektemp(i));
    end
end
temp_beta_1= temp_beta_1(~all(temp_beta_1 == 0, 2), :);
beta0=mean(temp_beta_1(:,2));
temp_up=zeros(307,2);
temp_down=zeros(307,2);
for i=1:303
    if weektemp(i)-weektemp(i+4)>0
        temp_down(i,1)=weektemp(i);
        temp_down(i,2)=iabeta(i);
    end
end
temp_down= temp_down(~all(temp_down == 0, 2), :);
for i=1:303
    if weektemp(i)-weektemp(i+4)<0
        temp_up(i,1)=weektemp(i);
        temp_up(i,2)=iabeta(i);
    end
end
temp_up= temp_up(~all(temp_up == 0, 2), :);


scatter(temp_up(:,1),temp_up(:,2)./ppval(pp,temp_up(:,1)), 'filled');
hold on
scatter(temp_down(:,1),temp_down(:,2)./ppval(pp,temp_down(:,1)), 'filled');
hold on
plot(-15:30,ones(46,1))
xlim([-8,26])
xlabel('temp（℃）')
ylabel("ratio")
% mean(temp_down(:,2)./beta0)
