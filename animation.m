close all

%data = solutionNonlinearLQR.y;
%Time = solutionNonlinearLQR.x; 

data = solutionNonlinearPP.y;
Time = solutionNonlinearPP.x; 


timeSeries = timeseries(data, Time);
dt = 0.05;
timeSeries = timeSeries.resample(timeSeries.Time(1):dt:timeSeries.Time(end));

rows_ = size(timeSeries.Time);
rows_ = rows_(1);

%
f1 = figure;
% animaci lze zrychlit pomoc� �pravy kroku, nap?. p?i 1:3:rows_ se
% vykresluje pouze ka�d� t?et� vzorek
for row = 1:1:rows_    

    alpha = timeSeries.Data(2, 1, row);
    Dalpha = timeSeries.Data(4, 1, row);

    xc = timeSeries.Data(1, 1, row);
    Dxc = timeSeries.Data(3, 1, row);
   
    %poloha kyvadla
    [xp, yp] = pol2cart(alpha-pi/2, L_p);
    %% vykreslovani
    cla
    hold on
    grid on
    
    axis equal
    %xlim([xc*1-0.1-L_p*0.5, xc*1+0.1+L_p*0.5])
    xlim([min(timeSeries.Data(1,1,:))-0.5, max(timeSeries.Data(1,1,:)+0.5)]);
    ylim([-L_p*1.5, L_p*1.5]);
    
    quiver( xc, 0,...
        xp, yp,...
        'Color', 'Black',...
        'LineWidth', 2,...
        'Marker', '*');

    time = timeSeries.Time(row);
    title("T = " + time + "   x = " + xc);

   drawnow
end
