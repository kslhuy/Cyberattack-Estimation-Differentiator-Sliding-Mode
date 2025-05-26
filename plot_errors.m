%% Results
figure(1)
set(gcf,'color','w');
subplot(3,1,1)
plot(t,ya11_h,'r-.',t,x(1,:),'b--','LineWidth',1.5);
% hold on
% plot(t,x_hat(1,:),'r-.',t,x(1,:),'b--','LineWidth',1.5);


grid on
xlabel('{Time} ${[sec]}$','Interpreter','LaTeX','Fontsize',12);
ylabel('$$\delta_i$$','Interpreter','LaTeX','Fontsize',12);
d1 = legend('$$\hat{\delta_i}$$','$$\delta_i$$');
d1.Interpreter = "latex";
d1.FontSize = 18;
d11 = title('Relative position between host and preceding vehicle and its estimation');
d11.Interpreter = "latex";
d11.FontSize = 10;


subplot(3,1,2)
set(gcf,'color','w');
% plot(t,ya12_h,'r-.',t,x(2,:),'b--','LineWidth',1.5);
plot(t,x_hat(2,:),'r',t,x(2,:),'b--','LineWidth',1.5);

grid on

xlabel('{Time} ${[sec]}$','Interpreter','LaTeX','Fontsize',12);
ylabel('\boldmath$\dot{\delta}_i$','Interpreter','LaTeX','Fontsize',12);
d2 = legend('$\dot{\hat{\delta}_i}$','$$\boldmath{\dot{\delta}_i}$$');
d2.Interpreter = "latex";
d2.FontSize = 18;
d12 = title('Relative velocity between host and preceding vehicle and its estimation');
d12.Interpreter = "latex";
d12.FontSize = 10;


subplot(3,1,3)
set(gcf,'color','w');
% plot(t,ya13_h,'r-.',t,x(3,:),'b--','LineWidth',1.5);
plot(t,x_hat(3,:),'r',t,x(3,:),'b--','LineWidth',1.5);

grid on
xlabel('{Time} ${[sec]}$','Interpreter','LaTeX','Fontsize',12);
ylabel('\boldmath$\ddot{\delta}_i$','Interpreter','LaTeX','Fontsize',12);
d3 = legend('$$\ddot{\hat{\delta}_i}$$','$$\boldmath{\ddot{\delta}_i}$$');
d3.Interpreter = "latex";
d3.FontSize = 18;
d13 = title('Relative acceleration between host and preceding vehicle and its estimation');
d13.Interpreter = "latex";
d13.FontSize = 10;


figure(2)
set(gcf,'color','w');
rgb1 = [0.4660 0.6740 0.1880];
rgb2 = [0 0.4470 0.7410];
% true
plot(t(2:end),nu,'Color',rgb1,'LineWidth',3);grid on
hold on;
% ouput one
plot(t, nu_h ,'Color',rgb2,'LineWidth',3);
% % ouput 2
plot(t, nu_h2 ,'LineWidth',3);
% 
moyen = (nu_h2 +  nu_h)/2;
plot(t,moyen ,'LineWidth',3);


grid on
xlabel('{Time} ${[sec]}$','Interpreter','LaTeX','Fontsize',12);
ylabel('$$f_{c}$$','Interpreter','LaTeX','Fontsize',12);
ylim([-30 30 ]);
d4 = legend('$$f_{c}$$' , '$\hat{f}_c1$','${\hat{f}_{c2}}$','$$f_{moy}$$');
% d4 = legend('$$f_{c}$$' , '$\hat{f}_c1$','${\hat{f}_{c2}}$','$$f_{moy}$$');



d4.Interpreter = "latex";
d4.FontSize = 18;
d14 = title('Cyber-attack and its estimation');
d14.Interpreter = "latex";
d14.FontSize = 10;



ax = axes;

set(ax,'units','normalized','position',[0.2,0.2,0.2,0.2])
box(ax,'on')
plot(t,nu_h,'Color',rgb1,'LineWidth',1.5);grid on
hold on;
plot(t,[nu nu_h(end)] ,'Color',rgb2,'LineWidth',1.5);grid on
set(ax,'xlim',[0,2],'ylim',[-85,20])