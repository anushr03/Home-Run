% Anush Rathod
% 01/12/2021
% ECE 202, Project 2, Phase 1
% Hitting a home run, with air resistance
% Phase 1 - Comparing the analytical solution to the numeric, without drag

clear; clf;

% ------- Given Information -------

v0mph = 112; % exit velocity in mph
phi0deg = 32;    % launch angle in degrees

x0 = 0; y0 = 0; % starting position coordinates for the ball

g = 10; % gravitational constant, in m/s

% ------- Setting up variables -------

mph2mps = 5280 * 12 * 2.54 / 100 / 3600;   % mph to m/s conversion
deg2rad = pi()/180;   % degrees to radians

v0 = v0mph * mph2mps;
phi0 = phi0deg * deg2rad;

v0x = v0*cos(phi0);   % x-component of v0
v0y = v0*sin(phi0);   % y-component of v0

% ----- compute some useful characteristics of trajectory -----

tH = v0y/g;    % time to reach max. height
tLand = 2*tH;   % time to land (time of flight)

H = tH * v0y/2;   % max. height
R = v0x * tLand ;  % range

m2ft = 3.3; % constant to convert from m to ft

R_ft = R*m2ft;    % ROUGH conversion from m to ft


% ----- set up a time array, compute x(t), y(t) analytically -----

tmin = 0; tmax = tLand; 
N = 2000;   % intervals

t = linspace(tmin, tmax, N+1);   % time array, connects x(t) with y(t)

% ------ Analytical Solution -------

xt = (x0 + v0x*t)*m2ft; % finding xt and converting to ft
yt = (y0 + v0y*t - (1/2)*g*t.^2)*m2ft;  % finding yt and converting to ft


% ------ Numerical Solution -------

m  = 0.145; % mass of baseball, in kg

dt = (tmax-tmin)/N;
y = zeros(1, N+1);   % initialize y(t)

x(1) = x0;
y(1) = y0;
vx = v0x;
vy = v0y;

for n = 1:N
    Fnet_x = 0;
    ax = Fnet_x/m;
    x(n+1) = x(n) + vx*dt + (1/2)*ax*dt^2;
    vx = vx + ax*dt;
    
    Fnet_y = -m*g;
    ay = Fnet_y/m;
    y(n+1) = y(n) + vy*dt + (1/2)*ay*dt^2; 
    vy = vy + ay*dt;
end

x_ft = x*m2ft; % converting x(t) to feet
y_ft = y*m2ft; % converting y(t) to feet

% ------ checks -------

checkSum_x = sum(abs(x_ft-xt))
checkSum_y = sum(abs(y_ft-yt))

% -------- plotting graph ---------

plot(xt, yt, x_ft, y_ft, 'LineWidth', 2)
ylim([0 140])
ax = gca; ax.FontSize = 15; ax.GridAlpha = 0.5;

xlabel('x (ft)', 'FontSize', 18)   % *** use ft for Project 2 ***
ylabel('y (ft)', 'FontSize', 18)

title({'ECE 202, Project 2, Phase 1, Trajectory of a baseball', ...
    'no drag, analytic vs. numeric solution'}, 'FontSize', 22)

legend({'analytic (behind numeric)', 'numeric'}, ...
    'FontSize', 18)

grid on

