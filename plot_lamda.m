%% Plot adaptive parameters lambda
figure;
set(gcf,'color','w');

subplot(2,1,1)
plot(t, lambda11, 'r', t, lambda12, 'g', t, lambda13, 'b', t, lambda14, 'k', 'LineWidth', 1.5);
legend('\lambda_{11}', '\lambda_{12}', '\lambda_{13}', '\lambda_{14}');
% xlabel('Time (s)');
ylabel('\lambda values');
title('Adaptive parameters for output 1');
grid on;

subplot(2,1,2)
plot(t, lambda21, 'r', t, lambda22, 'g', t, lambda23, 'b', 'LineWidth', 1.5);
legend('\lambda_{21}', '\lambda_{22}', '\lambda_{23}');
xlabel('Time (s)');
ylabel('\lambda values');
title('Adaptive parameters for output 2');
grid on;

% subplot(3,1,3)
% plot(t, z14, '--k', t, z23, '--b', 'LineWidth', 1.5);
% legend('z_{14}', 'z_{23}');
% xlabel('Time (s)');
% ylabel('z values');
% title('Auxiliary signals z_{14} and z_{23}');
% grid on;
