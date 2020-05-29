%%  Zadané hodnoty
clc
clear all 
close all
addpath('functions') % toto pøidá slo¾ku functions do prohlédávaných

p = getParameters();
s = tf('s');
% initializeModel();

%%  Model
    X_operating = [0 0 0 0]';
    [A, B] = AB(X_operating, 0);
    Co = [1 0 0 0; 0 0 1 0]; % observed outputs | measuring xc and alpha
    Cr = [1 0 0 0]; % reference outputs
    Cf = eye(4); % fully observed outputs
    D = [0; 0];
    Df = [0;0;0;0];

    P = ss(A,B,Co,D); %plant
    Po = ss(A,B,Cf,Df); %fully observed plant

%% Váhové filtry   
W1 = 1/makeweight(100,[10 1],0.1);
W2 = 0*s + 0.1;
W3 = 1/makeweight(0.1,[10 1],100);
%% Synteza H-inf kontroleru
[K,CL,gamma] = mixsyn(Po,W1,W2,W3);

loopsens_ = loopsens(Po,K);
S = loopsens_.Si;
T = loopsens_.Ti;

[svH, wH] = sigma(S);
Ms = max(svH(1,:));
[svH, wH] = sigma(T);
Mt = max(svH(1,:));
fprintf("Ms = %f dB \n Mt = %f dB\n",mag2db(Ms), mag2db(Mt));
%%
figure
p = sigmaplot(S,'b',T,'g',K*S,'r',gamma/W1, 'b-.',gamma/W2,'r-.', gamma/W3, 'g-.');
popt = getoptions(p);
popt.Title.String = 'Sensitivity functions';
popt.Grid = 'on';
setoptions(p,popt)
legend("S", "Gamma/W1", "T", "Gamma-inverted Weight 3");

%% Navrh estimatoru
    %observability check
    obsv_ = obsv(A,Co);
    obsvr = rank(obsv_);

    % process noise covariance matrix
    Vd = diag([5 10 5 10]);
    % noise covariance matrix
    Vn = [0.02 0; 0 0.05];
    
    L = lqe(A,Vd,Co,Vd,Vn);
    
    Cf = eye(4);
    KF = ss((A-L*Co), [B L], eye(4), 0);

%% Nastaveni pocatecnich hodnot
%pocatecni stav
X0 = [0,0,pi*2/16,0]'; %x, alpha, dx, dalpha
r = 0; %reference

simulationTime = 15;
dt = 0.005; %samplovaci perioda

%predalokace poli pro data
X = zeros(4, simulationTime/dt); %skutecny stav
X(:,1) = X0;
Xest = zeros(4,simulationTime/dt); %estimovany stav
Xest(:,1) = X0 - X_operating;
T = zeros(1,simulationTime/dt);   %spojity cas
U = zeros(1,simulationTime/dt);   %vstupy
U(1) = 0;
UX = zeros(length(K.A), simulationTime/dt); %vnitrni stav kontroleru
u = 0;
R = zeros(1,simulationTime/dt); %reference r
R(1) = r;
D = zeros(2,simulationTime/dt); %poruchy
Y = zeros(2,simulationTime/dt); %mereni
Y(:,1) = X0(1:2) - X_operating(1:2);

computingTimes = zeros(simulationTime/dt, 1);

d = [0 0];
d1T = 0;
d1t = 0;
d1a = 0;
d2T = 0;
d2t = 0;
d2a = 0;

%% Simulace
hbar = waitbar(0,'Simulation Progress');
tic
disp("1000 samples = " + 1000*dt + "s");
for k = 1:simulationTime/dt
    
    % Generovani pozadovane reference
%     if rand(1) > 0.999      
%         r = (2*rand(1)-1)*0.50;
%     end         

    %% Estimace stavu X
    y_est = Co * X(:,k);
    y_msr = Y(:,k);
    y_err = y_msr - y_est;
    Dxe = KF.A*Xest(:,k) + KF.B*[u ;y_msr];
%     disp("est: " + y_est + "  msr: " + y_msr + "  y_err: " + y_err)
    Xest(:,k+1) = Xest(:,k) + Dxe*dt; % Euler method
    %% Regulace
    w = [r 0 0 0]';
    e = -X(:,k);
    t = linspace(k*dt, (k+1)*dt);
    [utmp, t, UXtmp] = lsim(K,e'.*ones(length(t),4),t, UX(:,k));
    UX(:,k+1) = UXtmp(end);
%     Dux = K.A*UX(:,k) + K.B*e;
%     UX(:,k+1) = UX(:,k) + Dux*dt;
    u = utmp(end);
    u = min(12, max(-12, u));
   %% Simulace
    
    %"spojite" reseni v intervalu dt, uklada se pouze konecny stav 
    [ts, xs] = ode45(@(t, X_) pendCartC_d(X(:,k),u,d'), [(k-1)*dt k*dt], X(:,k));
    
    X(:,k+1) = xs(end,:)';
    T(k+1) = ts(end);
    R(k+1) = r;
    U(k+1) = u;
    
    % mereni Y
%         Y(:,k+1) = Co * xs(end,:)' + [ 0.02*sqrt(Vn(1,1))*randn(1), 0.003*sqrt(Vn(2,2))*randn(1) ]' ...
%                - [X_operating(1), X_operating(3)]';
    Y(:,k+1) = Co * xs(end,:)' - [X_operating(1), X_operating(3)]';

        

    %% Vizualizace
     
    %progress meter a vypocetni cas na 1000 vzorku
    if (mod(k,1000)==0) 
        disp("Computing time: " + toc)
        disp(k + "/" + simulationTime/dt);
        tic
    end

    waitbar(k*dt/simulationTime,hbar);   
end

close(hbar);

sol.X = X;
sol.R = R;
sol.Xest = Xest+X_operating;
sol.T = T;
sol.U = U;
sol.UX = UX;
sol.D = D;
sol.Y = Y;
%vytiskne øe¹ení
sol

save('results/ResultsHinf.mat', 'sol');
