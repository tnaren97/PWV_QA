baseDir = 'D:\PWV\PWV_data\Volunteers\Tarun_dQdA_Testing_110521\012_PWV_CartBH_AAo_PROSP_2VPS_81FRAMES_3p0ACC';
cd(baseDir);

cd('2DPC_QA_Analysis\AscAo');
load('FIGURE2_WORKSPACE.mat');

%% FLOW AREA WAVES
times = pcmrTimes(1:40);
smoothFlow = smoothdata(flow(1:40),'gaussian',4);
smoothArea = smoothdata(area(1:40),'gaussian',4);
IN = in(1:40);

flowFig = figure; flowFig.Position = [100 100 1500 800];
ax = gca;
ax.FontSize = 17; 

yyaxis left; plot(times,smoothFlow,'LineWidth',1.4);
hold on; scatter(times,smoothFlow,42,[0 0.4470 0.7410]);
hold on; scatter(times(IN),smoothFlow(IN),74,[0 0.4470 0.7410],'filled');
ylabel('Flow Rate (mL/s)','FontSize',24);
ylim([-60 360]);
yyaxis right; plot(times,smoothArea,'LineWidth',1.4);
hold on; scatter(times,smoothArea,42,[0.8500 0.3250 0.0980]);
hold on; scatter(times(IN),smoothArea(IN),74,[0.8500 0.3250 0.0980],'filled');
ylabel('Area (cm^2)','FontSize',24);
ylim([3.2 5.4]);
title('Ascending Aorta - Flow and Area Curves','FontSize',32); 
xlabel('Time (ms)','FontSize',24);
xlim([15 420]);


%% QA PLOT
x = smoothArea;
y = smoothFlow;
QAfig = figure; QAfig.Position = [100 100 1200 800];

scatter(x,y,72,times); 
ax = gca;
ax.FontSize = 20; 

c = colorbar;
c.Label.String = 'Time (ms)';
c.FontSize = 24;
title('Ascending Aorta - QA Plot','FontSize',32); 
xlabel('Area (cm^2)','FontSize',24);
ylabel('Flow (mL/s)','FontSize',24);
%free2 = drawfreehand;
%in = inpolygon(x,y,free2.Position(:,1),free2.Position(:,2));
%systolePts2 = [x(in); y(in)]';
%[coef,stats] = polyfit(systolePts2(:,1),systolePts2(:,2),1);

%Linear regression
%minArea = min(systolePts2(:,1));
%maxArea = max(systolePts2(:,1));
%xq = linspace(minArea-0.1,maxArea+0.1,100);
%yq = coef(1)*xq + coef(2); %y = mx+b
%delete(free2)
hold on; 
scatter(systolePts2(:,1),systolePts2(:,2),72,'k','x');
plot(xq,yq,'k','LineWidth',1.4); 
%str = ['Slope = ' num2str(coef(1)) '\newlineIntcpt = ' num2str(coef(2))];
%text(min(xq),max(yq)-0.1*max(yq),str);
hold off;
%DescAo_PWV = coef(1);
%disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);