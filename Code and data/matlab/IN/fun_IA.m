function fx=fun_IA(x)
global I_il
S0=x(1);
E0=x(2);
I0=x(3);
R0=x(4);
y0=[S0,E0,I0,R0,0];
tend=300*7;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
sol=ode45(@seir_IA,[0 tend-1],y0,options,x);
x1=linspace(0,tend-1,300); 
y1=deval(sol,x1);
weeknewcase=y1(5,2:end)-y1(5,1:end-1);
f1=weeknewcase-I_il(1:300);
fx=norm(f1);
end
