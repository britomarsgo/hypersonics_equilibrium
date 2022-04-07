% Plotting Graphics

Tc=[410.4966 664.8873 1.38E+03 3.00E+03 4.97E+03];
M=[4 7 12 19 25];
Ttn=[410.7834	646.6786	1.25E+03	2.43E+03	3.11E+03]

figure;
pos = get(gcf, 'Position');
plot(M,Tc,'k','linewidth',1.8);hold on
plot(M,Ttn,'red','linewidth',1.8);
ylabel('T2 (K)','FontSize',12);
xlabel('Mach Number','FontSize',12);
legend('Variação Temperatura Caloricamente Perfeito','Variação Temperatura em Equil.Químico','Location','Northwest')
grid on;
% axis([0.1 1 0.1 7*10^4 ])
set(gca,'YTick',[Tc(1),Tc(2),2000,Ttn(5),Tc(5),])
set(gca,'XTick',[4,7,25])
saveas(gcf,'T2versusmach.pdf')


pc=[6.11E+03	1.39E+04	3.60E+04	8.63E+04	1.47E+05];
M=[4 7 12 19 25];
ptn=[6.11E+03	1.38E+04	3.53E+04	8.31E+04	1.38E+05]

figure;
pos = get(gcf, 'Position');
plot(M,pc,'k','linewidth',1.8);hold on
plot(M,ptn,'red','linewidth',1.8);
ylabel('P2 (Pa)','FontSize',12);
xlabel('Mach Number','FontSize',12);
legend('Variação da Pressão Caloricamente Perfeito','Variação da Pressão em Equil.Químico','Location','Northwest')
grid on;
% axis([0.1 1 0.1 7*10^4 ])
set(gca,'YTick',[pc(1),pc(2),pc(5)])
set(gca,'XTick',[4,7,25])
saveas(gcf,'p2vsmach.pdf');

rc=[0.0519	0.0728	0.0909	0.1002	0.1033];
M=[4 7 12 19 25];
rtn=[0.0518	0.074	0.0975	0.1201	0.1461]

figure;
pos = get(gcf, 'Position');
plot(M,rc,'k','linewidth',1.8);hold on
plot(M,rtn,'red','linewidth',1.8);
ylabel('(\rho_{2}) (Kg/m³)','FontSize',12);
xlabel('Mach Number','FontSize',12);
legend('Variação da densidade Caloricamente Perfeito','Variação da densidade em Equil.Químico','Location','Northwest')
grid on;
% axis([0.1 1 0.1 7*10^4 ])
set(gca,'YTick',[rc(1),rc(2),rc(5),rtn(5)])
set(gca,'XTick',[4,7,25])
saveas(gcf,'rho2vsmach.pdf');


