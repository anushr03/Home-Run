% Anush Rathod
% 01/12/2021
% ECE 202, Project 2, Phase 5
% Hitting a home run, with air resistance
% Phase 1 - Comparing the analytical solution to the numeric, without drag
% phase 2 - Comparing drag and drag off trajectories
% phase 3 - exporting data to Excel and finding time of flight, maximum
%           height and range of the baseball
% phase 4 - Answering the same questions from phase 3, in MATLAB and
%           comparing results with Excel
% phase 5 - Exploring the result


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
R = v0x * tLand;  % range

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
C = input("Enter a value for drag coefficient C: ");
A = 0.00426; % cross sectional area, in m^2
            %(source-http://spiff.rit.edu/richmond/baseball/traj/traj.html)
p = 1.225; % air density, in kg/m^3 
          %(source - http://spiff.rit.edu/richmond/baseball/traj/traj.html)

dt = (tmax-tmin)/N;
y = zeros(1, N+1);   % initialize y(t)

x(1) = x0;
y(1) = y0;
vx = v0x;
vy = v0y;

export(1,:) = ["Time t (s)", "x (ft)","y (ft)"];
nRow = 2;
dt_temp = tmin; % temp variable to be used in matrix formation

for n = 1:N
    v = sqrt(vx^2 + vy^2);
    Fnet_x = 0 - (1/2)*C*p*A*v*vx;
    ax = Fnet_x/m;
    x(n+1) = x(n) + vx*dt + (1/2)*ax*dt^2;
    vx = vx + ax*dt;
    
    Fnet_y = -m*g - 0 - (1/2)*C*p*A*v*vy;
    ay = Fnet_y/m;
    y(n+1) = y(n) + vy*dt + (1/2)*ay*dt^2; 
    vy = vy + ay*dt;
    
    %  matrix formation to export to excel 
    
    export(nRow, :) = [dt_temp, (x(n)*m2ft), (y(n)*m2ft)];
    dt_temp = dt_temp + dt;
    nRow = nRow + 1;
    
    % finding time of flight and range
    
    if y(n)/y(n+1) < 0
        time = t(n); % getting the corresponding value from time array
        range = x(n); % getting the corresponding value from x array
        vF = v; % getting the velocity at which ball lands
    end
    
end

x_ft = x*m2ft; % converting x(t) to feet
y_ft = y*m2ft; % converting y(t) to feet

% -------- answers to various questions --------

time_of_flight = time % in sec
maximum_height_ft = max(y_ft) % in ft
range_ft = range*m2ft % in ft

% The answers check out with the answers found in the excel sheet from
% phase 3

% ----- final velocity of the ball ------

Velocity_final = vF/mph2mps % in miles/hr
Energy_lost_J = (1/2)*m*(vF^2 - v0^2) 

% ------ checks -------

checkSum_x = sum(abs(x_ft-xt))
checkSum_y = sum(abs(y_ft-yt))

% -------- plotting graph ---------

plot(xt, yt, x_ft, y_ft, 'LineWidth', 2)
ylim([0 140])
ax = gca; ax.FontSize = 15; ax.GridAlpha = 0.4; ax.MinorGridAlpha=0.5;
xlabel('x (ft)', 'FontSize', 18) 
ylabel('y (ft)', 'FontSize', 18)

title({'ECE 202, Project 2, Phase 5, Trajectory of a baseball', ...
    'drag vs no drag'}, 'FontSize', 22)
str1 = sprintf("drag, C = %g", C);
legend('no drag', str1, 'FontSize', 18)

grid on
grid minor

% ----- export data to a file -----

writematrix(export, 'trajectory.csv')

% ------- Calculating percent errors to the given diagram -------

Percent_error_timeFlight = (time - 5.7)/5.7
Percent_error_maxHeight = (maximum_height_ft - 114)/114
Percent_error_Range = (range_ft - 446)/446

% c) As we can see with the percent error, the calculated values are
% pretty close to the given values and they compare pretty well with them. 



