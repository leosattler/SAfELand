function classfigsmaker(map, outputFolder, fileName)
%==========================================================================
% Function to generate plot (.png) for final susceptibility map.
%
% Input types: (array, string, string). 
% map = map of values (for watershed or scars, either q/t ratio or shalstab classes)
% outputFolder = output folder location
% fileName = plot file name (suffix)
%==========================================================================
% Creating array for plot (map clone)
arr=map;
%-------------------------------------------------------------------------
% Retrieving pixels inside watershed 
bacia=map>-9999;
%--------------------------------------------------------------------------
% Retrieving pixels from uncond. stable and UNstable classes
classesIncond=map==10 | map==-10;
%--------------------------------------------------------------------------
% Retrieving pixels that ARE NOT from unconditional classes
classesInterm=bacia & ~classesIncond;  
%--------------------------------------------------------------------------
% Retrieving pixels from each intermediary classes
class3=map>=-9.9 & map<=-3.1 & classesInterm; % -9.9 <= log q/T <= -3.1 
class4=map>-3.1 & map<=-2.8 & classesInterm;  % -3.1 <  log q/T <= -2.8 
class5=map>-2.8 & map<=-2.5 & classesInterm;  % -2.8 <  log q/T <= -2.5
class6=map>-2.5 & map<=-2.2 & classesInterm;  % -2.5 <  log q/T <= -2.2
class7=map>-2.2 & map<=+9.9 & classesInterm;   % -2.2 < log q/T <= +9.9
%--------------------------------------------------------------------------
% Chancing pixel values for each corresponding class
arr(map==10)=2;
arr(class7)=3;
arr(class6)=4;
arr(class5)=5;
arr(class4)=6;
arr(class3)=7;
arr(map==-10)=8;
%--------------------------------------------------------------------------
% Counting amount of pixels inside each class for percentage display in the
% plot
total=sum((sum(arr>-9999)));
c2=num2str(round(100*sum(sum(arr==2))/total,2));
c3=num2str(round(100*sum(sum(arr==3))/total,2));
c4=num2str(round(100*sum(sum(arr==4))/total,2));
c5=num2str(round(100*sum(sum(arr==5))/total,2));
c6=num2str(round(100*sum(sum(arr==6))/total,2));
c7=num2str(round(100*sum(sum(arr==7))/total,2));
c8=num2str(round(100*sum(sum(arr==8))/total,2));
% Finalizing string with percentage infos
c2=strcat(c2,'%');
c3=strcat(c3,'%');
c4=strcat(c4,'%');
c5=strcat(c5,'%');
c6=strcat(c6,'%');
c7=strcat(c7,'%');
c8=strcat(c8,'%');
%--------------------------------------------------------------------------
% Declaring figure 
fig=figure('Visible', 'off', 'DefaultAxesPosition', ...
    [0.004, 0.025, 0.875, 0.95]);
%--------------------------------------------------------------------------
% Configuring figure size
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperPosition', [0 0 26 13]); %x_width=26cm y_width=13cm
%--------------------------------------------------------------------------
% Altering pixel values outside scar to 0
%arr(map==-9999)=0;
%--------------------------------------------------------------------------
% Plotting suceptibility map array
image(arr);
%--------------------------------------------------------------------------
% Declaring colorbar
cmap=[1 1 1; 0 0 1; 0 0.5 1; 0 1. 1.; 0.5 1. 0.5; 1. 1. 0; 1. 0.6 0; 1 0 0];
colormap(cmap);
cb = colorbar; 
cb.LineWidth=0.1;
cb.FontSize=11.2;
% Configuring colorbar labels - text and position
set(cb,'xtick', [1.5 [2.4 2.6] [3.4 3.6] [4.4 4.6] [5.4 5.6] [6.4 6.6] [7.4 7.6] [8.4 8.6]], 'xticklabel', ... 
    {'Sem Classe', ...
    sprintf(strcat(c2,'% \n','Incond. Estavel')) ...
    sprintf(strcat(c3,'% \n','log(q/T) > -2.2')), ...
    sprintf(strcat(c4,'% \n','-2.5 < log(q/T) < -2.2')), ...
    sprintf(strcat(c5,'% \n','-2.8 < log(q/T) < -2.5')), ...
    sprintf(strcat(c6,'% \n','-3.1 < log(q/T) < -2.8')), ... 
    sprintf(strcat(c7,'% \n','log(q/T) < -3.1')), ...
    sprintf(strcat(c8,'% \n','Incond. Instavel')) ...
    }, 'ylim', [2 9]);
%--------------------------------------------------------------------------
% Shutting down X and Y axis in the plot
axis off;
%--------------------------------------------------------------------------
% Saving plot 
print(fig,[outputFolder '\' fileName '.png'],'-dpng','-r400');
print(fig,[outputFolder '\' fileName '.jpeg'],'-djpeg','-r400');
print(fig,[outputFolder '\' fileName '.pdf'],'-dpdf','-r400');
end