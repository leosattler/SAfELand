function [nscars, AllScars] = scarsidentifier(map)
%==========================================================================
% Function to find most representative value of log(q/T) inside scar and
% apply it uniformly to the entire scar polygon. The routine scans the grid
% in search for pixels that are inside scars by checking jumps from NoData
% pixel to one with a value (inside scar). It attributes an id to all pixels 
% pertaining to a scar each time a scar is found. If a pixel with no value 
% is found and it is neighboring one already valued in a previous scan, 
% the id value is copied. At the end, we use the  functions uniformscarmin 
% to assign an id to possibly unvalued pixels, since a scar shape may be 
% so that the performed scans might still miss a pixel. 
%
% Input types: (array)
% map = map of values (for watershed or scars, either q/t ratio or shalstab classes)
%==========================================================================
% Finding number of cols and rows in grid
[rows, columns]=size(map);
%--------------------------------------------------------------------------
% Creating uniform grid
AllScars=map;
AllScars(map~=-9999)=0;
%--------------------------------------------------------------------------
% Initiating loop (n=1 is id of first found scar)
n=1;
for i=2:rows-1
    for j=2:columns-1
        if AllScars(i,j) - AllScars(i,j-1) == 9999 % Pixel turns from NoData to Data: scar was found! (id=n) 
            if AllScars(i-1,j) ~= -9999 || ...   % Case in which neighboring pixels were already identified (id = n of row above)
               AllScars(i-1,j-1) ~= -9999 || ... 
               AllScars(i-1,j+1) ~= -9999 || AllScars(i,j-1) ~= -9999
                AllScars(i,j) = max([AllScars(i-1,j) AllScars(i-1,j-1) AllScars(i-1,j+1) AllScars(i,j-1)]);
            else                                 % Case in which scar has no id, attributing id (id= n +1, to count new scar)
                AllScars(i,j) = n;               
            end
        elseif AllScars(i,j) - AllScars(i,j-1) > -9999 && AllScars(1,j) - AllScars(i,j-1) < 0 % Still inside same scar: keep repeating id 
            AllScars(i,j) = AllScars(i,j-1);                                                  
        elseif AllScars(i,j) - AllScars(i,j-1) <= -9999 % Scan exited scar and returned to NoData pixel: next id will be n+1
            n = n + 1;                                  
        end
    end
end
%--------------------------------------------------------------------------
% Due to the shape of some scars (specially alongated ones), some remain
% with more than one id for is pixels. To make it uniform the minimum n of any 
% neighboring pixel is attributed to the entire scar with the function unofrmscarsmin.
UnifMapFinal=uniformscarsmin(AllScars);
% Finding all indexes (length(listofscars) = the number of scars inside
% watershed!)
listofscars=intersect(UnifMapFinal(UnifMapFinal~=-9999),UnifMapFinal(UnifMapFinal~=-9999));
nscars=length(listofscars);
AllScars=map; % recreating grid for final classification
% Finally, we attribute the mode of vale log(q/T) inside each scar
% Reads: UnifMap(UnifMap=cicatriz) = mode(grid(UnifMap=cicatriz))
for i=1:length(listofscars)
    AllScars(UnifMapFinal==listofscars(i))=i;
end
end