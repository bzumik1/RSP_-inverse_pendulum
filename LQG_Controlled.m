%%  Zadan� hodnoty
% clc
clear all 
close all
addpath('functions') % toto p�id� slo�ku functions do prohl�d�van�ch

%vytvo�en� struct pro v�echny zadan� hodnoty syst�mu

p = getParameters();
initializeModel();
%%  Model
    X_operating = [0 0 pi 0]';
    [A, B] = AB(X_operating, 0);
    Co = [1 0 0 0; 0 0 1 0]; % observed outputs | measuring xc and alpha
    Cr = [1 0 0 0]; % reference outputs
    D = [0; 0];
    
    % Adding an error integrator into the state space description
    % The controller's objective is to follow the reference r
    % with the variable x_c, where x_c is the position of the cart
    % The controlled plant already has an integrator, the controller's
    % integrator purpose is to reject constant disturbances, not to
    % remove steady-state error. 
    Ah = [A, zeros(length(A),1);
         -Cr, 0];
    Bh = [B; 0];
    Cho = [1 0 0 0 0; 0 0 1 0 0];
    Dh = [0;0];
    
    Gr = ss(Ah,Bh,Cho,Dh);
    Ge = ss(A,B,Co,D);
%%  Navrh regulatoru
    ctrb_ = ctrb(Ah,Bh);
    ctrbr = rank(ctrb_);

    Q = diag([1 0.1 1 1 10]);
    R = 0.1;
    [K,S,e] = lqr(Ah,Bh,Q,R);

%% Navrh estimatoru
    obsv_ = obsv(A,Co);
    obsvr = rank(obsv_);

    % process noise covariance matrix
    Vd = diag([5 10 5 10]);
    % noise covariance matrix
    Vn = [0.02 0; 0 0.05];
    
%     [kest,L,P] = kalman(Ge,Vd,Vn);
     L = lqe(A,Vd,Co,Vd,Vn);
%   L = (lqr(A', Co', Vd, Vn))';

    Cf = eye(4);
    KF = ss((A-L*Co), [B L], eye(4), 0);
   
%% Nastaveni pocatecnich hodnot
%pocatecni stav
X0 = [0,0,pi,0]'; %x, alpha, dx, dalpha
%pozadovany stav x_cr
r = 0.5;

%nastaveni solveru
options = odeset();

simulationTime = 2.5e1;
dt = 0.01; %samplovaci perioda
kRefreshPlot = 10; %vykresluje se pouze po kazdych 'kRefreshPlot" samplech
kRefreshAnim = 2; % ^

%predalokace poli pro data
X = zeros(4, simulationTime/dt); %skutecny stav
X(:,1) = X0;
DKsi = zeros(1,simulationTime/dt); %odchylka od reference
Ksi = zeros(1,simulationTime/dt);
Xest = zeros(4,simulationTime/dt); %estimovany stav
Xest(:,1) = X0 - X_operating;
Ts = zeros(1,simulationTime/dt);   %spojity cas
U = zeros(1,simulationTime/dt);   %vstupy
U(1) = 0;
u = 0;
R = zeros(1,simulationTime/dt); %reference r
R(1) = r;
D = zeros(2,simulationTime/dt); %poruchy
Y = zeros(2,simulationTime/dt); %mereni
Y(:,1) = X0(1:2) - X_operating(1:2);

%pocatecni porucha
d = [0 0]'; %vektor poruch
d1T = 0; %celkove trvani poruchy 1
d1t = 0; %jak dlouho je porucha 1 momentalne aktivni
d1a = 0; %amplituda poruchy 1

%% Simulace
hbar = waitbar(0,'Simulation Progress');
tic
disp("1000 samples = " + 1000*dt + "s");
for k = 1:simulationTime/dt
    
    % Generovani pozadovane reference
    if rand(1) > 0.999      
        r = (2*rand(1)-1)*0.50
    end               

    %% Regulace
    %definice vstupu a saturace do <-6,6>
    
      e = [ Xest(:,k) ; Ksi(k) ]; %vektor v rozsirenem stavovem prostoru
%       e(1) = -r + e(1); % signs are reversed because they get flipped by -K_lqr
      u = -K * e;
      u = min(12, max(-12, u));
     
      if mod(k,500)==0
          for i = 1:length(e)
              str = sprintf("k%d * e%d = %f", i,i,-K(i)*e(i));
              disp(str);
          end
        str = sprintf("u = %f \n",u);
        disp(str);
      end
    %% Estimace stavu X

    y_est = Co * X(:,k);
    y_msr = Y(:,k);
    y_err = y_msr - y_est;
    Dxe = KF.A*Xest(:,k) + KF.B*[u ;y_msr];
%     disp("est: " + y_est + "  msr: " + y_msr + "  y_err: " + y_err)
    Xest(:,k+1) = Xest(:,k) + Dxe*dt; % Euler method
    %% Simulace
    
    %"spojite" reseni v intervalu dt, uklada se pouze konecny stav 
    [ts, xs] = ode45(@(t, X_) pendCartC_d(X(:,k),u,d), [(k-1)*dt k*dt], X(:,k), options);
    
    X(:,k+1) = xs(end,:)';
    Ts(k+1) = ts(end);
    U(k+1) = u;
    R(k+1) = r;
    D(:,k+1) = d;
    
    % mereni Y
    Y(:,k+1) = Co * xs(end,:)' + ...
               - [X_operating(1), X_operating(3)]';

    Dksi(k+1) = r - Y(1,k+1);
    Ksi(k+1) = Ksi(k) + dt*Dksi(k);
        

    %% Vizualizace
     
    %refresh plotu
    if(mod(k+1,kRefreshPlot)==1)
%         plotRefresh(Ts,X,Xest+X_operating,R,U,D,Y,k,kRefreshPlot);
    end
       
    %progress meter a vypocetni cas na 1000 vzorku
    if (mod(k,1000)==0) 
        disp("Computing time: " + toc)
        disp(k + "/" + simulationTime/dt);
        tic
    end

    waitbar(k*dt/simulationTime,hbar);
    
    %%
    if abs(X(1, k)) > 10 || abs(X(3,k)-pi) > pi
        break;
    end
    
end
close(hbar);



sol.X = X;
sol.Xest = Xest+X_operating;
sol.T = Ts;
sol.U = U;
sol.R = R;
sol.D = D;
sol.Y = Y;
sol.dt = dt;

sol

save('results/ResultsLQG.mat', 'sol');
