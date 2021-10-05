%% Housekeeping
clc
clear 
close all

%% About the Project:
%Goals:
    % Use ode45 to integrate a projectile trajectory
    % Find the position trajectory for a three-phase water rocket
    % Adjust parameters to achieve 80 +/- 1 m landing
    
%Steps:
    % 1) Set constants
    % 2) Create a function that describes the three phases
    % 3) Pass it to ode45 and integrate
    % 4) Plot the trajectory and thrust
    
%Assumptions:
    % Compressible and isentropic properties of flow
    % No other external forces such as wind
    % Constant gravity
    % 2D motion
    
%Units: SI

%% Constants
Pa      = 83426.56;  % Ambient pressure [Pa]
rhoAir  = 0.961;     % Density of ambient air [kg/m^3]
rhoH2O  = 1000;      % Density of water [kg/m^3]

gamma   = 1.4;       % ratio of specific heats for air
R       = 287;       % Air gas constant [J / kg / k]
g       = 9.81;      % gravitational constant
Cd      = 0.8;       % discharge coefficient

massB        = 0.15;      % Bottle empty mass
volumeB      = 0.002;     % Bottly empty volume
TempAir_i    = 300;       % Air initial temp [K]
Dthroat      = 0.021;     % Throat diameter  [m]
DB           =  0.105;    % Bottle diameter  [m]

% Calculate areas:
A_t = pi / 4 * Dthroat^2;
A_b = pi / 4 * DB^2;

% Test Parameters:
theta        = pi / 4 ;      % initial angle of rocket
CD           = 0.5;       % coefficient of drag
volumeH2O_i  = 0.001;     % Water initial volume within bottle [m^3]
Pabs         = 428164;  % Initial air absolute pressure in bottle (gage + atm) [Pa]

%Initial Conditions:
vx0         = 0;
vz0         = 0;
x0          = 0;
z0          = 0.25;
volumeAir_i = volumeB - volumeH2O_i; % Air initial volume
mAir0       = (Pabs * volumeAir_i) / (R * TempAir_i); % Inital air mass
m0          = massB + mAir0 + volumeH2O_i * rhoH2O; % Rocket inital mass
length      = 0.5;    % test stand length
tspan       = [0 5];   % Integration time
%% ODE45:
Y0 = [x0 z0 vx0 vz0 m0  mAir0 volumeAir_i]'; %inital state
% Constants vector:
C  = [Pabs, Pa, rhoAir, rhoH2O, gamma, R, g, CD, Cd, massB, volumeH2O_i,...
      volumeB, TempAir_i, A_t, A_b, volumeAir_i, theta, length, mAir0];

options = odeset('Events', @eventfunction); % Ode constraint (check function)
[t, Y]  = ode45(@(t, Y) rocketEquation(t, Y, C), tspan, Y0, options);

x       = Y(:,1);
z       = Y(:,2);  
vx      = Y(:,3);
vz      = Y(:,4);
mass    = Y(:,5);
Airmass = Y(:, 6);
Airvol  = Y(:, 7);

%% Thrust Equation:
[imax, ~] = size(Airvol);
Pend      = Pabs * (volumeAir_i / volumeB)^gamma;
u = [1 1 1];
for i = 1:imax
   if (Airvol(i) < volumeB)
    Pressure(i,1) = Pabs * (volumeAir_i / Airvol(i))^gamma;
    Fthrust(i,1)  = 2 * Cd * A_t * (Pressure(i) - Pa); 
    u(1) = i;
   end
   
   if (Airvol(i) >= volumeB && Pend > Pa)
    Pressure(i,1) = Pend * (Airmass(i) / mAir0)^gamma;
    rhoA = Airmass(i) / volumeB;
    Tair = Pressure(i) / (rhoA * R);
    Pcritical(i,1) = Pressure(i) * (2 / (gamma + 1))^(gamma/(gamma-1));
    
      if (Pcritical(i) > Pa)
         Te    = (2 * Tair) / (gamma + 1);  % Exhaust temperature
         rho_e = Pcritical(i) / (R * Te);   % Exhaust air density
         Ve    = sqrt(gamma * R * Te); 
         Pe    = Pcritical(i);  
      elseif (Pcritical(i) < Pa)
         Mach  = sqrt(2 * (Pressure(i) / Pa)^((gamma-1)/gamma) - 2) / sqrt(gamma -1);
         Te    = Tair / (1 + (gamma-1)/2 * Mach^2);
         rho_e = Pa / (R * Te);
         Ve    = Mach * sqrt(gamma * R * Te); 
         Pe    = Pa;       % Exhaust pressure 
       end
     u(2) = i;  
    Fthrust(i) = Cd * rho_e * A_t * Ve^2 + A_t * (Pa - Pe);
   end
   
   if(Airvol(i) >= volumeB && Pressure(i) <= Pa)
       Fthrust(i) = 0;
       u(3) = i; 
   end
end
%% Plotting:

[maxHeight, imax] = max(z);
phase2i   = find(Fthrust == 0,1);  % Index when phase 3 starts
maxRange  = max(x);

% subplot(2,3,1) 
figure(1)
plot(x, z); hold on
caption1 = sprintf('z = %.2f m', maxHeight);
plot(x(imax), maxHeight, '.r'); text(x(imax) -17 , maxHeight, caption1), hold on
caption2 = sprintf('x = %.2f m', maxRange);
text(maxRange-5, 3, caption2)
xline(x(u(1)), 'g'); xline(x(phase2i), 'r')
xlabel(" Range [ m ]");   ylabel(" Height [ m ]");
legend("Projectile", "Max Height", "Phase 2", "Phase 3", 'Location', 'South');
title(" Projectile Motion of a Water Rocket");
grid on
hold on

figure(2)
% subplot(2,3,2)
plot(t(1:70), Fthrust(1:70));
xlabel(" Time [s] "); ylabel(" Thrust [N] "); hold on
xline(t(u(1)), 'g'); xline(t(phase2i), 'r')
legend("Thrust", "Phase 2", "Phase 3")
title (" Thrust envelope of the Water Rocket");
grid

figure(3)
plot(t(1:60), Airvol(1:60));
xlabel(" Time [s]");  ylabel("Volume of Air in Bottle Rocket [m^3]");
title("Bottle Rocket Air Volume Evolution");
grid

figure(4)
plot(t, vz);
xlabel("Time [s]");   ylabel("Z Velocity [m/s]");
title("Bottle Rocket Vertical Velocity Profile");
grid

figure(5)
plot(t, vx);
xlabel("Time [s]");   ylabel("X Velocity [m/s]");
title("Bottle Rocket Horizontal Velocity Profile");
grid

%% Functions:

function [dydt] = rocketEquation(t, Y, C) 

% Define State Vector:
x           = Y(1);           % Horizontal distance [m]
z           = Y(2);           % Vertical distance   [m]
vx          = Y(3);          
vz          = Y(4);
mass        = Y(5);           % Rocket mass         [kg]
massAir     = Y(6);
AirVolume   = Y(7);

% Redefine constants from C:
Pairi  = C(1);     Pa     = C(2);
rhoAmb = C(3);     rhoH2O = C(4);
gamma  = C(5);     R      = C(6);
g      = C(7);     CD     = C(8);
Cd     = C(9);     massB  = C(10);
vH2o_i = C(11);    volB   = C(12);
T0     = C(13);    A_t    = C(14); 
A_b    = C(15);    vAir_i = C(16);
theta  = C(17);    l      = C(18);
mAir_i = C(19); 


% Phase 1: Water Thrust

    % Final state of phase 1 (independent of time):
    Pend     = Pairi * (vAir_i / volB)^gamma;   %Pressure
    Tend     = T0 * (vAir_i / volB)^(gamma-1);  %Temperature
    
    
if(AirVolume < volB)
    % Volumetric rate of change:
    vdotAir = Cd * A_t * sqrt(2/rhoH2O * (Pairi * (vAir_i / AirVolume)^gamma - Pa));
    % Pressure change:
    P       = Pairi * (vAir_i / AirVolume)^gamma; 
    % Rocket mass flow rate:
    mdot    = -Cd * A_t * sqrt(2 * rhoH2O * (P - Pa)); 
    % Thrust:
    thrust  = 2 * Cd * A_t * (P - Pa); 
    mdotAir = 0; 
  
end
% Phase 2: Air thrust 

if (AirVolume >= volB && Pend > Pa)
    % Air pressure:
    P         = Pend * abs((massAir / mAir_i))^gamma; 
    % Air density:
    rhoAir    = massAir / volB;  
    % Air temperature:
    Tair      = P / (rhoAir * R);
    
    % Defining critical pressure:
    Pcritical = P * (2 / (gamma + 1))^(gamma/(gamma-1));

    if(Pcritical > Pa)
        Te    = (2 * Tair) / (gamma + 1);  % Exhaust temperature
        rho_e = Pcritical / (R * Te);   % Exhaust air density
        Ve    = sqrt(gamma * R * Te); disp('chocked')
        Pe    = Pcritical;              % Exhaust pressure
        
    elseif(Pcritical < Pa)
        Mach  = sqrt(2 * (P / Pa)^((gamma-1)/gamma) - 2) / sqrt(gamma -1);
        Te    = Tair / (1 + (gamma-1)/2 * Mach^2); 
        rho_e = Pa / (R * Te);disp('unchocked')
        Ve    = Mach * sqrt(gamma * R * Te); 
        Pe    = Pa;       % Exhaust pressure 
        
    end
    
    mdotAir  = -Cd * rho_e * A_t * Ve;   %%%% negative or positve
    mdot     =  mdotAir; 
    vdotAir  = 0;
    thrust   = -mdotAir * Ve + A_t * (Pa - Pe); 
end

% Phase 3: Free Fall 
if(AirVolume >= volB && P <= Pa)
    thrust  = 0; 
    mdot    = 0; 
    mdotAir = 0;
    vdotAir = 0;

end

% Rocket heading direction: z < length*sin(theta)
if(x <= l* cos(theta) && z<= 0.25 + (l*sin(theta)))
    hx = cos(theta);
    hz = sin(theta);
else
    hx = vx / sqrt(vx^2 + vz^2);
    hz = vz / sqrt(vx^2 + vz^2); 
end

% Drag equation:
Drag         = 0.5 * rhoAmb * (vx^2 + vz^2) * CD * A_b;
DragX        = Drag * hx; % Drag projection onto X
DragZ        = Drag * hz; % Drag projection onto Z

% Thrust Componenets:
thrustX      = thrust * hx;
thrustZ      = thrust * hz; 

% Acceleration:
% acceleration = (thrust - Drag + g * mass) / mass;
accX         = (thrustX - DragX) / mass;
accZ         = (thrustZ - DragZ - g * mass) / mass;


dydt(1) = vx;           % Velcity in the x direction 
dydt(2) = vz;           % Velocity in the z direction
dydt(3) = accX;         % Acceleration in th x direction
dydt(4) = accZ;         % Acceleration in the z direction
dydt(5) = mdot;         % Mass flow rate of rocket
dydt(6) = mdotAir;      % Mass flow rate of air
dydt(7) = vdotAir;      % Volumetric flow rate of air

dydt = dydt'; 
end

% Terminate ode45 when z = 0:
 function [position, isterminal, direction] = eventfunction(~, Y)
       position   = Y(2);  % Second element in state vector which is z
       isterminal = 1;
       direction  = 0;
        
 end

