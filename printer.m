function printer(indx, toadd)            
%==========================================================================
% Auxiliar function for visualizaton of simulation evolution and to save
% generated outpus
%
% Input types: (double, double). 
% indx = index representing the stage of parameter combination
% toadd = value to add to the internal function counter
%==========================================================================
% Index i is tha value to be altered in loop
i = indx;
% Case in which this is the last output; code i=-9999 is used to
% print the following message
if i==-9999
    pause(.2) 
    for t=1:3 % blink three times in different intervals
        clc;
        pause(.02*t*t)    
        fprintf('Finished! %s', ['|' [repmat('=',1,50) '='] '| ' num2str(100) '%']);
        fprintf('\n');     
        fprintf('Saving outputs...  100%%');
        fprintf('\n'); 
        fprintf('Finishing...');
        fprintf('\n');     
        pause(.2*t*t)
    end
else % Case in which counter is ongoing
    % Animated wheel has 9 stages. If counter passes, return to its initial 
    % state (if i=18, returns to 1, if i=28, returns to 2...)
    if i>9 
        i = i - 9*fix(i/9) + 1;
    end
    % Percentage representing progess (29 is max value arbitrarily defined,
    % that approaches number of operations during save of outputs)
    % Updates occur only in some cases 
    if i == 5
        percent = [num2str(round(100*(indx-1)/(29+toadd),2)) '%%'];
    elseif i == 6
        percent = [num2str(round(100*(indx-2)/(29+toadd),2)) '%%'];
    elseif i == 7
        percent = [num2str(round(100*(indx-3)/(29+toadd),2)) '%%'];
    elseif i == 8
        percent = [num2str(round(100*(indx-4)/(29+toadd),2)) '%%'];
    else
        percent = [num2str(round(100*indx/(29+toadd),2)) '%%'];
    end
    % Initiating loop for animation
    if i>=1 && i <=3     % updating progress
        s = '.  ';
	elseif i>=4 && i<=6
        s = '.. ';
	elseif i>=7 && i<=9 
        s = '...';
    end
    if 2 <= i && i <= 8
        clc;
        fprintf('Finished! %s', ['|' [repmat('=',1,50) '='] '| ' num2str(100) '%']);
        fprintf('\n');     
        fprintf(['Saving outputs' s ' ' num2str(percent)]);
        fprintf('\n');
    end
end
end  