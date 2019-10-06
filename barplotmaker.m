function barplotmaker(ClassBasin, ClassScars, outputFolder, fileName)
%==========================================================================
% Function to generate plot (.png) of shalstab classes' histogram with
% a priori probability curve.
%
% Input types: (array, array, string, string).
% ClassBasin = shalstab class map for watershed
% ClassScars = shalstab class map for scars
% outputFolder = output folder location
% fileName = plot file name (suffix)
%==========================================================================
% Retrieving % of pixels inside each shalstab class for watershed
c1=ClassBasin(ClassBasin<=-10 & ClassBasin<-9.9 & ClassBasin~=-9999);
c2=ClassBasin(ClassBasin>=-9.9 & ClassBasin<=-3.1 & ClassBasin~=-9999);
c3=ClassBasin(ClassBasin>-3.1 & ClassBasin<=-2.8 & ClassBasin~=-9999);
c4=ClassBasin(ClassBasin>-2.8 & ClassBasin<=-2.5 & ClassBasin~=-9999);
c5=ClassBasin(ClassBasin>-2.5 & ClassBasin<=-2.2 & ClassBasin~=-9999);
c6=ClassBasin(ClassBasin>-2.2 & ClassBasin<=9.9 & ClassBasin~=-9999);
c7=ClassBasin(ClassBasin>9.9 & ClassBasin>=10 & ClassBasin~=-9999);
% Creating list of percentages
BarBasin=[length(c1) length(c2) length(c3) length(c4) length(c5) length(c6) length(c7)];
totalBasin=length(ClassBasin(ClassBasin~=-9999));
%--------------------------------------------------------------------------
% Retrieving % of pixels inside each shalstab class for scars
cScars1=ClassScars(ClassScars<=-10 & ClassScars<-9.9 & ClassScars~=-9999);
cScars2=ClassScars(ClassScars>=-9.9 & ClassScars<=-3.1 & ClassScars~=-9999);
cScars3=ClassScars(ClassScars>-3.1 & ClassScars<=-2.8 & ClassScars~=-9999);
cScars4=ClassScars(ClassScars>-2.8 & ClassScars<=-2.5 & ClassScars~=-9999);
cScars5=ClassScars(ClassScars>-2.5 & ClassScars<=-2.2 & ClassScars~=-9999);
cScars6=ClassScars(ClassScars>-2.2 & ClassScars<=9.9 & ClassScars~=-9999);
cScars7=ClassScars(ClassScars>9.9 & ClassScars>=10 & ClassScars~=-9999);
% Creating list of percentages
BarScars=[length(cScars1) length(cScars2) length(cScars3) length(cScars4) length(cScars5) ...
    length(cScars6) length(cScars7)];
totalScars=length(ClassScars(ClassScars~=-9999));
%--------------------------------------------------------------------------
% Generating plots
fig3=figure('Visible', 'off', 'DefaultAxesPosition', ...
    [0.12, 0.26, 0.81, 0.73]);
set(fig3, 'PaperUnits', 'centimeters');
set(fig3, 'PaperPosition', [0 0 17 12]); %x_width=26cm y_width=13cm
% Plot with 2 y axis
y1=[100*(BarBasin/totalBasin)' 100*(BarScars/totalScars)']; % Barplot
x=[.72 2.17 3.57 4.97 6.37 7.77 9.17]; % Valores do eixo x sao arbitrarios (precisam ser consecutivos e conter 7)
y2=1*(totalBasin/totalScars)*(BarScars./BarBasin); % Plot estilo stem e linha
% Plotting
[hAx,~,plot2] = plotyy(x, y1, x, y2, 'bar', 'plot');
plot2.MarkerSize=6;
plot2.Marker='+';
% Plot configuration
ylabel(hAx(1),'Pixels [%]','FontSize', 14) % label eixo y da esquerda
ylabel(hAx(2),'%Scars/%Basin','FontSize', 14) % label eixo y da direita
plot2.LineWidth = 1.2; % line width
ylim(hAx(1),[0 100]);
% Generating legend (only for current plot)
legend('Basin', 'Scars','% Ratio','Location','north','Orientation','horizontal');
hold on
% Axis limits
ylim(hAx(1),[0 (max(max(y1))+10)]);
ylim(hAx(2),[0 (max(y2)+5)]);
% Adding texts
for i=1:7
    % Probability string info
    % Configuring text vertical position for probability curve (blue line) 
    % which axis lies on the right
    ax2texty=(((max(max(y1)))/(max(y2)+2))*y2(i));
    % Checking if text has overlap with other graph infos, e.g. 
    % bars and its text values (strings of percentage values)   
    if (ax2texty <= (y1(i,1)+5) && (y1(i,2)-5) <= ax2texty) || ...
            (ax2texty <= (y1(i,2)+5) && (y1(i,1)-5) <= ax2texty) && i > 1
        ax2texty = max(y1(i,:)) + 7;
    end
    text((x(i)-.05),ax2texty,num2str(round(y2(i),1)),'color',[0 0.45 0.74],'FontSize', 11);
    % Texts for watershed bars 
    if length([num2str(round(y1(i,1),1)) '%']) == 3
        text((x(i)-.52),(y1(i,1)+3.1),[num2str(round(y1(i,1),1)) '%'],'FontSize', 11);
    elseif length([num2str(round(y1(i,1),1)) '%']) == 4
        text((x(i)-.65),(y1(i,1)+3.1),[num2str(round(y1(i,1),1)) '%'],'FontSize', 11);
    elseif length([num2str(round(y1(i,1),1)) '%']) == 5
        text((x(i)-.74),(y1(i,1)+3.1),[num2str(round(y1(i,1),1)) '%'],'FontSize', 11);
    else
        text((x(i)-.41),(y1(i,1)+3.1),[num2str(round(y1(i,1),1)) '%'],'FontSize', 11);
    end
	% Texts for scars' bars
    text((x(i)+.07),(y1(i,2)+3.1),[num2str(round(y1(i,2),1)) '%'],'FontSize', 11);
end
% Plot configuration
plot2.LineWidth = 1; % line width
% x axis label configuration
set(gca,'xTick',x);
set(gca,'xticklabel',{'Uncond. Unstable'; 'log(q/T) < -3.1'; '-3.1 < log(q/T) < -2.8'; ... 
    '-2.8 < log(q/T) < -2.5'; '-2.5 < log(q/T) < -2.2'; 'log(q/T) > -2.2'; 'Uncond. Stable'}, ...
    'XTickLabelRotation',36,'FontSize', 13);
%--------------------------------------------------------------------------
% Saving plot
print(fig3,[outputFolder '\' 'histogram_' fileName '.png'],'-dpng','-r400');
print(fig3,[outputFolder '\' 'histogram_' fileName '.jpeg'],'-djpeg','-r400');
print(fig3,[outputFolder '\' 'histogram_' fileName '.pdf'],'-dpdf','-r400');
end