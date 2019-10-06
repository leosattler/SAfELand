function UnifMap = uniformscarsmin(map)
%==========================================================================
% Function to find minimum value of log(q/T) of scar and apply it uniformly
% to entire scar polygon. The script scans each pixel (i,j) to find minimum
% possible value within the 9 analyzed pixels: 8 neighboring pixels + 1
% (i,j):
%              
%                (i-1,j)
% (i-1,j-1) +       +       + (i-1,j+1)
%           
%           
%   (i-1,j) +       +       + (i,j+1)
%                 (i,j)
%
% (i-1,j+1) +       +       + (i+1,j+1)
%                (i,j+1)
%
% Input types: (array).
%
% map = map of values (for watershed or scars, either q/t ratio or shalstab classes)
% i = row index
% j = collumn index
%==========================================================================
% Finding number of rows and columns
[rows, columns]=size(map);
%--------------------------------------------------------------------------
% Creating an uniform grid
UnifMap = map;
%--------------------------------------------------------------------------
% Defyning -9999 values (NoData) as +9999 for minimum value check
UnifMap(UnifMap==-9999) = 9999;
%--------------------------------------------------------------------------
% Initiating loop (from left to right, from top to bottom)
for i=1:rows
    for j=1:columns
        % Pixel grid limits
        Li=i-1;
        Lf=i+1;
        Ci=j-1;
        Cf=j+1;
        % Avoiding pixels outside grid boundary:
        %(i or j = 0, i or j > size(grid))
        if Li==0       % on first row
            Li=1;
        end
        if Ci==0       % on first col
            Ci=1;
        end
        if Lf>rows     % on last  row
            Lf=rows;
        end
        if Cf>columns  % on last col
            Cf=columns;
        end
        % Finding smallest value (avoinding NoData)
        if UnifMap(i,j)<9999
            UnifMap(i,j)=min(min(UnifMap(Li:Lf,Ci:Cf)));
        end
    end
end
%--------------------------------------------------------------------------
% Initiating loop (from left to right, from bottom to top)
for i=-1*(-rows:-1)
    for j=1:columns
        % Pixel grid limits
        Li=i-1;
        Lf=i+1;
        Ci=j-1;
        Cf=j+1;
        % Avoiding pixels outside grid boundary:
        %(i or j = 0, i or j > size(grid))
        if Li==0       % on first row
            Li=1;
        end
        if Ci==0       % on first col
            Ci=1;
        end
        if Lf>rows     % on last  row
            Lf=rows;
        end
        if Cf>columns  % on last col
            Cf=columns;
        end
        % Finding smallest value (avoinding NoData)
        if UnifMap(i,j)<9999
            UnifMap(i,j)=min(min(UnifMap(Li:Lf,Ci:Cf)));
        end
    end
end
%--------------------------------------------------------------------------
% Initiating loop (from right to left, from top to bottom)
for i=1:rows
    for j=-1*(-columns:-1)
        % Pixel grid limits
        Li=i-1;
        Lf=i+1;
        Ci=j-1;
        Cf=j+1;
        % Avoiding pixels outside grid boundary:
        %(i or j = 0, i or j > size(grid))
        if Li==0       % on first row
            Li=1;
        end
        if Ci==0       % on first col
            Ci=1;
        end
        if Lf>rows     % on last  row
            Lf=rows;
        end
        if Cf>columns  % on last col
            Cf=columns;
        end
        % Finding smallest value (avoinding NoData)
        if UnifMap(i,j)<9999
            UnifMap(i,j)=min(min(UnifMap(Li:Lf,Ci:Cf)));
        end
    end
end
%--------------------------------------------------------------------------
% Initiating loop (from right to left, from bottom to top)
for i=-1*(-rows:-1)
    for j=-1*(-columns:-1)
        % Pixel grid limits
        Li=i-1;
        Lf=i+1;
        Ci=j-1;
        Cf=j+1;
        % Avoiding pixels outside grid boundary:
        %(i or j = 0, i or j > size(grid))
        if Li==0       % on first row
            Li=1;
        end
        if Ci==0       % on first col
            Ci=1;
        end
        if Lf>rows     % on last  row
            Lf=rows;
        end
        if Cf>columns  % on last col
            Cf=columns;
        end
        % Finding smallest value (avoinding NoData)
        if UnifMap(i,j)<9999
            UnifMap(i,j)=min(min(UnifMap(Li:Lf,Ci:Cf)));
        end
    end
end
%--------------------------------------------------------------------------
% Turning NoData pixels into -9999 again
UnifMap(UnifMap==9999)=-9999;
end
