yyaxis left
% tempf = linspace(1,40,40);
% tempf(1:15) = AscAo_flow(1:15);
% tempf(16:35) = 0;
% tempf(36:40) = AscAo_flow(16:20);
% plot(circshift(AscAo_flow,15))
selected = circshift(asystole,-7);
flows = circshift(AscAo_flow, -7);
plot(pcmr_times,flows)
% plot(pcmr_times(1:18),flows(1:18))

hold on;
plot(pcmr_times(selected), flows(selected), 'k*')
ylabel('Flow (cm^3/s)')
% tempa = linspace(1,40,40);
% tempa(1:15) = AscAo_area(1:15);
% tempa(16:35) = 0;
% tempa(36:40) = AscAo_area(16:20);
yyaxis right
% plot(circshift(AscAo_area,15))
plot(pcmr_times,circshift(smoothdata(AscAo_area_int),-7))
% plot(circshift(smoothdata(AscAo_area, 'gaussian', 5),15))
ylabel('Area (cm^2)')
title('Ascending Aorta - Flow and Area curves')
xlim([1 860])
xlabel('Time (ms)')
saveas(gcf,'AscAo_AreaFlow_plot');
% figure; scatter(AscAo_area, AscAo_flow, [], linspace(1, length(AscAo_area), length(AscAo_area)))
% colorbar()
% 
% figure
% yyaxis left
% plot(circshift(AscAo_flow_new,15))
% ylabel('Flow')
% yyaxis right
% plot(circshift(AscAo_area_new,15))
% ylabel('Area')