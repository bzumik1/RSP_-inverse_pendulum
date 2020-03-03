clc
close all
addpath('functions')
addpath('gif')
addpath('Results')


data = load('Results1.mat');
data = data.sol;

samples = length(data.X);

Xs = data.X;
Xest = data.Xest;
Ts = data.T;
U = data.U;
Wx = data.Wx;
D = data.D;
Y = data.Y;

if (isfield(data, 'computingTimes'))
    computingTimes = data.computingTimes;
    figure("Name", "Computing times")
    bar(Ts(1:end-1), computingTimes);
end


bonked_k = data.bonked_k

kRefreshPlot = 10; %vykresluje se pouze po kazdych 'kRefreshPlot" samplech
kRefreshAnim = 5; % ^

        %animRefresh(Ts,Xs,[],1);
        %gif('NLMPC_Swingup_dt10_sawtooth_highV.gif')
for k = 2:1:samples-1
    %% Vizualizace
    if(mod(k,kRefreshPlot)==0)
       %plotRefresh(Ts,Xs,[],[],U,D,Y,k,kRefreshPlot);
    end
    
    if(mod(k,kRefreshAnim)==0)
        animRefresh(Ts,Xs,Wx,k);
        title(k)
        %gif
    end
        
    if(mod(k,10000)==0)
        disp("Time for 10000 samples:" + toc)
        tic
    end
    

end
