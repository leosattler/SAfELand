function [outputArea,X,Y,Xtmp,Ytmp]=integral(mapBasin, mapScars, resolution)
%==========================================================================
% Auxiliar function for the suceptibility map validation based on the
% greatest integral value of X vs Y curve, where:
% areaValue = integral 
% Xtmp = n# of pixels of each class inside watershed
% Ytmp = n# of pixels of each class inside scars
% X = sum of pixels in each class for watershed (in %)
% Y = sum of pixels in each class for scars (in %)
% They are equivalent to (example with 10 pixels and 3 classes, with 2 
% pixels in class 1; 3 pixels in class 2; 5 pixels in class 3):
% Xtmp = [2 5 3]
% X = [2 5 10]/10 == [2 2+3 2+3+5]/10 == [0.2 0.5 1]
%
% Input types: (array, array, double). 
% mapBasin = map of values for watershed
% mapScars = map of values for scars
% resolution = Number of digits after decimal point to be truncated when 
% note: the greater the resolution, the greater is the number of defined 
% classes and, consequently, greater the number of points for validation
% curve. A well defined (and behaved) curve is ideal for better determining
% its integral and, for such, more reliable is the validation process.
%==========================================================================
% Truncating DEM class values based on chosen resolution
mapBasin=round(mapBasin,resolution);
mapScars=round(mapScars,resolution);
%--------------------------------------------------------------------------
% Creating list with all existing classes 
classes=intersect(mapBasin(mapBasin~=-9999),mapBasin(mapBasin~=-9999)); 
% Defining n# of classes based on given resolution 
nclass=length(classes);
% Creating temporary list (tmp) and final list for X axis
Xtmp=ones(1,nclass); 
X=zeros(1,length(Xtmp));
% Creating temporary list (tmp) and final list for Y axis
Ytmp=zeros(1,nclass);
Y=zeros(1,length(Xtmp));
% Making list with all different values of NoData inside watershed and
% scars
mB=sort(mapBasin(mapBasin~=-9999));
mS=sort(mapScars(mapScars~=-9999));
%--------------------------------------------------------------------------
% Counting classes inside watershed
c=1; % Class counter
C=0; % Global class counter
j=1; % Class index
for i=2:length(mB)
    if mB(i-1) - mB(i) == 0  % If current pixel class is same as previous,
        c = c + 1;           % counter is increased for next check (c=c+1) and
        Xtmp(j) = c;         % its value is used for given class inside Xtmp list.
                            
    else                     % If current pixel class differs from previous, 
        C = C+Xtmp(j);       % global counter C is increased for X list,
        X(j)=C;              % its value is used for given class inside X list,
        c = 1;               % counter c is reset (c=1) and  
        j = j+1;             % index is increased for next class (j=j+1).
    end
end
% Finalizando lista de contagens acumuladas das classes na bacia - X
X(end)=sum(Xtmp(1:end));                     
X=X/X(end);
%--------------------------------------------------------------------------
% Counting classes inside scars 
classesScars=intersect(mS,mS); % List of classes inside scars
Ytmp0=ones(1,length(classesScars)); % Creating temporary list (tmp)(0).
                                    % It will be used to generate main Ytmp 
                                    % list ahead.
c=1; % Class counter
j=1; % Class index
for i=2:length(mS)
    if mS(i-1) - mS(i) == 0  % If current pixel class is same as previous,
        c = c + 1;           % counter is increased for next check (c=c+1) and
        Ytmp0(j) = c;        % its value is used for given class inside Ytmp0  list.
                            
    else                     % If current pixel class differs from previous, 
        c = 1;               % counter is reset (c=1) and
        j = j+1;             % index is increased for next class (j=j+1).
    end
end
% Creating final list for Ytmp based on Ytmp0, which has info of amount
% of classes inside watershed, but zeros were not considered (classes inside
% watershed but not on scars). Here this is corrected. First a list is
% created for positions of classes in scars in the list of all classes 
% inside watershed (classes inside scars are a subset of watershed classes):
[~,i1,~]=intersect(classes,classesScars); 
for i=1:length(i1) 
    indx=i1(i);    
    Ytmp(indx)=Ytmp0(i); 
end                      
% Creating list of the accumulated count of classes inside scars - Y list
Cy=0;
for i=1:length(Xtmp)
    % Y (ONLY SCARS):
    Cy=Cy+Ytmp(i);
    Y(i)=Cy;
end
Y=Y/Y(end);    
%--------------------------------------------------------------------------
% Applying trapz function to calculate integral 
outputArea=trapz(X, Y);
end