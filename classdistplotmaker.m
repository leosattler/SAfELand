function classdistplotmaker(THETA1, FLOWACC1, ClassBasin, outputFolder, fileName)
%==========================================================================
% Function to generate plot (.png) for class distribution of watershed
% pixels as a function of the log of accumulation area (y axis) and slope
% (x axis) according to shalstab susceptibility map.
%
% Input types: (array, array, array, string, string).
% thetamap1 = file name for watershed slopes [rad]
% flowaccmap1 = file name for watershed contributing area [m2] 
% ClassBasin = watershed shalstab classes map
% outputFolder = output folder location
% fileName = plot file name (suffix)
%==========================================================================
% Colors to be sued for shalstab classes classification
color_red = [1 0 0];  % red
color_orange = [1. 0.5 0];  % orange
color_yellow = [1.0000 1.0000 0];  % yellow
color_green = [0.5000 1.0000 0.5000];  % green
color_lblue = [0 1. 1];  % light blue
color_blue = [0 0.5 1];  % blue
color_dblue = [0 0 1];  % dark blue
%--------------------------------------------------------------------------
% Defining (column) vector of colors
colors = [color_red; color_orange; color_yellow; color_green; color_lblue; color_blue; color_dblue];
%--------------------------------------------------------------------------
% Defining maps for x and y and dot size (for plot)
x=THETA1; % (x>0 indicates dots inside watershed)
y=log10(FLOWACC1); % (y axis = LOG(FLOWACC))
tamanhoponto=1; % dot size desired for graph
%--------------------------------------------------------------------------
% x
[x1,ind1]=sort(x(ClassBasin<=-10 & ClassBasin<-9.9 & x>0));
[x2,ind2]=sort(x(ClassBasin>=-9.9 & ClassBasin<=-3.1 & x>0));
[x3,ind3]=sort(x(ClassBasin>-3.1 & ClassBasin<=-2.8 & x>0));
[x4,ind4]=sort(x(ClassBasin>-2.8 & ClassBasin<=-2.5 & x>0));
[x5,ind5]=sort(x(ClassBasin>-2.5 & ClassBasin<=-2.2 & x>0));
[x6,ind6]=sort(x(ClassBasin>-2.2 & ClassBasin<=9.9 & x>0));
[x7,ind7]=sort(x(ClassBasin>9.9 & ClassBasin>=10 & x>0));
%--------------------------------------------------------------------------
% y
y1=y(ClassBasin<=-10 & ClassBasin<-9.9 & x>0);y1=y1(ind1);
y2=y(ClassBasin>=-9.9 & ClassBasin<=-3.1 & x>0);y2=y2(ind2);
y3=y(ClassBasin>-3.1 & ClassBasin<=-2.8 & x>0);y3=y3(ind3);
y4=y(ClassBasin>-2.8 & ClassBasin<=-2.5 & x>0);y4=y4(ind4);
y5=y(ClassBasin>-2.5 & ClassBasin<=-2.2 & x>0);y5=y5(ind5);
y6=y(ClassBasin>-2.2 & ClassBasin<=9.9 & x>0);y6=y6(ind6);
y7=y(ClassBasin>9.9 & ClassBasin>=10 & x>0);y7=y7(ind7);
%--------------------------------------------------------------------------
% dot color (this is redundant; however, it avoids erros related to the 
% existence of classes to be represented)
c1=ClassBasin(ClassBasin<=-10 & ClassBasin<-9.9 & x>0);
c2=ClassBasin(ClassBasin>=-9.9 & ClassBasin<=-3.1 & x>0);
c3=ClassBasin(ClassBasin>-3.1 & ClassBasin<=-2.8 & x>0);
c4=ClassBasin(ClassBasin>-2.8 & ClassBasin<=-2.5 & x>0);
c5=ClassBasin(ClassBasin>-2.5 & ClassBasin<=-2.2 & x>0);
c6=ClassBasin(ClassBasin>-2.2 & ClassBasin<=9.9 & x>0);
c7=ClassBasin(ClassBasin>9.9 & ClassBasin>=10 & x>0);
%--------------------------------------------------------------------------
% Reducing number of dots to the mean at each N values. This is to
% make a less dot-crowded graph.
N = 2;
% 1
X1=[];
Y1=[];
% Defining general indexes
i=0;
f=0;
for n=1:length(x1)/N
    i=1+f;
    f=n*N;
    X1(end+1)=mean(x1(i:f));
    Y1(end+1)=mean(y1(i:f));
end
%
% 2
X2=[];
Y2=[];
% Defining general indexes
i=0;
f=0;
for n=1:length(x2)/N
    i=1+f;
    f=n*N;
    X2(end+1)=mean(x2(i:f));
    Y2(end+1)=mean(y2(i:f));
end
%
% 3
X3=[];
Y3=[];
% Defining general indexes
i=0;
f=0;
for n=1:length(x3)/N
    i=1+f;
    f=n*N;
    X3(end+1)=mean(x3(i:f));
    Y3(end+1)=mean(y3(i:f));
end
%
% 4
X4=[];
Y4=[];
% Defining general indexes
i=0;
f=0;
for n=1:length(x4)/N
    i=1+f;
    f=n*N;
    X4(end+1)=mean(x4(i:f));
    Y4(end+1)=mean(y4(i:f));
end
%
% 5
X5=[];
Y5=[];
% Defining general indexes
i=0;
f=0;
for n=1:length(x5)/N
    i=1+f;
    f=n*N;
    X5(end+1)=mean(x5(i:f));
    Y5(end+1)=mean(y5(i:f));
end
%
% 6
X6=[];
Y6=[];
% Defining general indexes
i=0;
f=0;
for n=1:length(x6)/N
    i=1+f;
    f=n*N;
    X6(end+1)=mean(x6(i:f));
    Y6(end+1)=mean(y6(i:f));
end
%
% 7
X7=[];
Y7=[];
% Defining general indexes
i=0;
f=0;
for n=1:length(x7)/N
    i=1+f;
    f=n*N;
    X7(end+1)=mean(x7(i:f));
    Y7(end+1)=mean(y7(i:f));
end
%--------------------------------------------------------------------------
% Finding maximum and minimum values to be displayed in the plot (this
% will help plot configuration)
ymin=0.8*min(y(x>0));
ymax=1.2*max(y(x>0));
xmin=0.8*min(x(x>0));
xmax=1.2*max(x(x>0));
%--------------------------------------------------------------------------
% Artificially correcting possible class dots that cross other classes
% threshold value
X1 = X1( X1==X1 & X1>mean(x2) & X1>mean(x3) & X1>mean(x4) & X1>mean(x5) ...
    & X1>mean(x6) );
Y1 = Y1( X1==X1 & X1>mean(x2) & X1>mean(x3) & X1>mean(x4) & X1>mean(x5) ...
    & X1>mean(x6) );
X7 = X7( X7==X7 & X7<mean(X3) & X7<mean(X3) & X7<mean(X4) & X7<mean(X5) ...
    & X7<mean(X6) );
Y7 = Y7( X7==X7 & X7<mean(X3) & X7<mean(X3) & X7<mean(X4) & X7<mean(X5) ...
    & X7<mean(X6) );
%--------------------------------------------------------------------------
% Initiating plot
fig1=figure('Visible', 'off', 'DefaultAxesPosition', ...
    [0.075, 0.145, 0.91, 0.83]);
set(fig1, 'PaperUnits', 'centimeters');
set(fig1, 'PaperPosition', [0 0 14 9]); %x_width=26cm y_width=13cm
%--------------------------------------------------------------------------
% Checking if shalstab class is present in the watershed and checking
% its legend for display
legend1=0;
legend2=0;
legend3=0;
legend4=0;
legend5=0;
legend6=0;
legend7=0;
%--------------------------------------------------------------------------
% Plotting only those dots that are present in the watershed
if ~isempty(c1)
    plot(X1,Y1, '.', 'MarkerSize', tamanhoponto, 'Color', colors(1,:)); hold on;
    legend1=1;
end
if ~isempty(c2)
    plot(X2,Y2, '.', 'MarkerSize', tamanhoponto, 'Color', colors(2,:)); hold on;
    legend2=1;
end
if ~isempty(c3)
    plot(X3,Y3, '.', 'MarkerSize', tamanhoponto, 'Color', colors(3,:)); hold on;
    legend3=1;
end
if ~isempty(c4)
    plot(X4,Y4, '.', 'MarkerSize', tamanhoponto, 'Color', colors(4,:)); hold on;
    legend4=1;
end
if ~isempty(c5)
    plot(X5,Y5, '.', 'MarkerSize', tamanhoponto, 'Color', colors(5,:)); hold on;
    legend5=1;
end
if ~isempty(c6)
    plot(X6,Y6, '.', 'MarkerSize', tamanhoponto, 'Color', colors(6,:)); hold on;
    legend6=1;
end
if ~isempty(c7)
    plot(X7,Y7, '.', 'MarkerSize', tamanhoponto, 'Color', colors(7,:)); hold on;
    legend7=1;
end
%--------------------------------------------------------------------------
% Drawing vertical lines and text 'x=...'
if ~isempty(X1) && ~isempty(X7) % Checando se classes incondicionais estao presentes na bacia
    % First line
    % Method to find 1st line position, avoiding class gaps or 
    % minimum/maximum dot values that cross other class thresholds:
    % 1- find minimum value of all intermediary classes (leftmost values)
    % 2- find maximum value among these (rightmost value from step 1)
    % 3- find maximum value of uncond. stable class (rightmost class value)
    % 4- find minimum value (leftmost) among values from steps 2 and 3
    pos1 = min([max([min(X2) min(X3) min(X5) min(X6)]) max(X7)]);
    % Testing if line position is still not ideal    
    if max([min(X2) min(X3) min(X5) min(X6)]) < mean(X7)
        pos1 = max(X7);
    end
    line([pos1 pos1], [ymin ymax],'Color','k'); hold off;
    t1 = double(1.01*pos1);
    t2 = double(.9*ymax);
    text(t1, t2, strcat('x1=',num2str(pos1)), 'FontSize', 12);
    % Second line
    % Method with same logic of 1st line but for opposite graph extreme:    
    % 1- find maximum value of all intermediary classes (rightmost values)
    % 2- find minimum value among these (leftmost value from step 1)
    % 3- find minimum value of uncond. UNstable class (leftmost class value)
    % 4- find maximum value (rightmost) among values from steps 2 and 3
    pos2 = max([min(X1) min([max(X2) max(X3) max(X5) max(X6)])]); 
    % Testing if line position is still not ideal  
    if min([max(X2) max(X3) max(X5) max(X6)]) > mean(X1)
        pos2 = min(X1);
    end
    line([pos2 pos2], [ymin ymax],'Color','k'); hold on;
    t1 = double(1.01*pos2);
    t2 = double(.8*ymax);
    text(t1,t2,strcat('x2=',num2str(pos2)), 'FontSize', 12);
end
%--------------------------------------------------------------------------
% Defining x and y axis limits
axis([xmin xmax ymin ymax]);
%ax = gca;
%ax.FontSize = 13; 
%--------------------------------------------------------------------------
% Defining x and y axis labels
xlabel('Slope [rad]','FontSize', 13);
ylabel('Log(acc)','FontSize', 13);
%--------------------------------------------------------------------------
% Checking ideal legend for the analized parameters
legendstrings=cell(1,(legend1+legend2+legend3+legend4+legend5+legend6+legend7));
for i=1:(legend1+legend2+legend3+legend4+legend5+legend6+legend7)
    if legend1==1  % if legend=1, legend will be added
        legendstrings(i)={strcat('-10  — -9.9 (', num2str(100*length(c1)/length(x(x>0)),2), '%)')};
        legend1=0; % changing to 0 so it won't be added again
    elseif legend2==1
        legendstrings(i)={strcat('-9.9 — -3.1 (', num2str(100*length(c2)/length(x(x>0)),2), '%)')};
        legend2=0;
    elseif legend3==1
        legendstrings(i)={strcat('-3.1 — -2.8 (', num2str(100*length(c3)/length(x(x>0)),2), '%)')};
        legend3=0;
    elseif legend4==1
        legendstrings(i)={strcat('-2.8 — -2.5 (', num2str(100*length(c4)/length(x(x>0)),2), '%)')};
        legend4=0;
    elseif legend5==1
        legendstrings(i)={strcat('-2.5 — -2.2 (', num2str(100*length(c5)/length(x(x>0)),2), '%)')};
        legend5=0;
    elseif legend6==1
        legendstrings(i)={strcat('-2.2 —  9.9 (', num2str(100*length(c6)/length(x(x>0)),2), '%)')};
        legend6=0;
    elseif legend7==1
        legendstrings(i)={strcat(' 9.9 —  10  (', num2str(100*length(c7)/length(x(x>0)),2), '%)')};
        legend7=0;        
    end
end
%--------------------------------------------------------------------------
[~,icons]=legend(legendstrings); % Calling legend
for i=8:length(icons)
    icons(i).MarkerSize=22;
end
%lgd.FontSize = 14;
%--------------------------------------------------------------------------
% Saving plot
print(fig1,[outputFolder '\' 'class_distribution_' fileName '.png'],'-dpng','-r400');
print(fig1,[outputFolder '\' 'class_distribution_' fileName '.jpeg'],'-djpeg','-r400');
print(fig1,[outputFolder '\' 'class_distribution_' fileName '.pdf'],'-dpdf','-r400');
end