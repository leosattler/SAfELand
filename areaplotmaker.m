function areaplotmaker(QTBasin, QTScars, resolution, outputFolder, fileName)
%==========================================================================
% Function to generate plot (.png) of curve (+ integral value) for the
% area validation method.
%
% Inputs types: (array, array, double, string, string).
% QTBasin = Watershed's q/t ratio
% QTScars = Scars' q/t das cicatrizes
% resolution = Number of digits after decimal point to be truncated when 
% storing pixel's q/t value in the q/t class intervals
% outputFolder = output folder location
% fileName = plot file name (suffix)
%==========================================================================
% Area directly from integral function
[areaValue,X,Y,~,~]=integral(QTBasin, QTScars, resolution);
%--------------------------------------------------------------------------
% Generating figure
fig2=figure('Visible', 'off', 'DefaultAxesPosition', ...
    [0.1, 0.13, 0.86, 0.8]);
set(fig2, 'PaperUnits', 'centimeters');
set(fig2, 'PaperPosition', [0 0 15 12]); 
plot(X,Y,'b+');
hold on;
plot(X,Y,'k-');
legend(['area = ' num2str(areaValue)],'Location','east');
xlabel('Cumulative % of Classes - Basin','FontSize', 13);
ylabel('Cumulative % of Classes - Scars','FontSize', 13);
%--------------------------------------------------------------------------
% Saving figure
print(fig2,[outputFolder '\' 'integral_' fileName '.png'],'-dpng','-r400');
print(fig2,[outputFolder '\' 'integral_' fileName '.jpeg'],'-djpeg','-r400');
print(fig2,[outputFolder '\' 'integral_' fileName '.pdf'],'-dpdf','-r400');
end