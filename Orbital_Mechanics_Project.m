clear, close, clc;

%% Assignment
R = [8228, 389, 6888]; % Geocentric position [km]
V = [-0.7, 6.6, -0.6]; % Geocentric velocity [km/s]

radius_E = 6378; % Radius of Earth [km]
radius_J = 71490; % Radius of Jupiter [km] 

mu_sun = 1.327E+11; % Gravitational parameter for ther Sun [km^3/s^2]
mu_jupiter = 126.686E+06; % Gravitational parameter for Jupiter [km^3/s^2]
mu_earth = 398600; % Gravtiational paramter for Earth [km^3/s^2]

R_E = 149.6E+06; % Distance between the Sun and Earth 
R_J = 778.6E+06; % Distance between the Sun and Jupiter [km]

t_E = 1; % Period of Earths orbit [Earth years]
t_J = 11.86; % Period of Jupiters orbit [Earth years]

I_sp_solid = 290; % Specific impulse for solid propellant [sec]
I_sp_liquidH = 455; % Specific impulse for liquid hydrogen propellant [sec]
g_0 = 9.807E-03; % Sea-level standard acceleration of gravity [km/s^2]

%% Inital Calculations
R_mag = norm(R); % Magnitude of geocentric position [km]
V_mag = norm(V); % Magnitude of geocentric velocity [km/s]
v_r = dot(R, V); % Radial vecloity [km/s]
V_E = sqrt(mu_sun/R_E); % Helocentric velocity of Earth
V_J = sqrt(mu_sun/R_J); % Helocentric velocity of Jupiter

%% Part A - Orbital Elements
h = cross(R, V); % Specific angular momentum vector [km^2/s]
h_mag = norm(h); % Magnitude of specific angular momentum [km^2/s]
C = cross(V, h) - mu_earth*(R/R_mag); % Laplace vector [km^3/s^2]
e = C/mu_earth; %Eccentricity vector [unitless]
e_mag = norm(e); % Magnitude of eccentricity [unitless]
i = acosd(h(1, 3)/h_mag); % Inclination angle [degrees]
N = [-1*h(1, 2), -1*-1*h(1, 1), 0]; % Node line vector [unitless] note: cross product won't work?? -> [-1*h(1, 2),-(-h(1, 1)), 0]
N_mag = norm(N); % Magnitude of Node line [unitless];
Omega = [acosd(N(1, 1)/N_mag), 360 - (acosd(N(1, 1)/N_mag))]; % Right ascension of the ascending node [degrees]
omega = [acosd(dot(N, e)/(N_mag*e_mag)), 360 - acosd(dot(N, e)/(N_mag*e_mag))]; % Argument of perigee [degrees]
theta = [acosd(dot(e, R)/(e_mag*R_mag)), 360 - acosd(dot(e, R)/(e_mag*R_mag))]; % True anomaly [degrees]

% Displaying Orbital Elements
fprintf('Part A: Orbital elements\n') % Printing introductory sentence
% Specific anglular momentum
fprintf('h = %.2f km^2/s\n', h_mag); % Printing the value for specific angular momentum

% Eccentricity
fprintf('e = %.5f\n', e_mag); % Printing the value for eccentricity

% Inclination angle
if i > 0 && i <= 90 % Determining if incination indicates a prograde orbit
    type = 'Prograde';
    angle_1 = 0;
    angle_2 = 90;
elseif i > 90 && i <= 180
    type = 'Retrograde'; % Determining if inclination indicates a retrograde orbit
    angle_1 = 90;
    angle_2 = 180;
end
fprintf('i = %.2f%c  //  %d%c < i \x2264 %d%c --> %s\n', i, char(176), angle_1, char(176), angle_2, char(176), type); % Printing the value for inclination angle and the reasoning

% Right ascension of the ascending node
if N(1, 2) > 0 % Determining of RAAN based on the value of the y-comp of Node line
    type = 'N_y > 0';
    angle_1 = 0;
    angle_2 = 180;
    if Omega(1, 1) >= 0 && Omega (1, 1) < 180 % Determining which value in the matrix is the smaller one
        angle = Omega(1, 1);
    else 
        angle = Omega(1, 2);
    end 
else 
    type = 'N_y < 0'; % Determining RAAN based on the value of y-comp of Node line
    angle_1 = 180;
    angle_2 = 360;
    if Omega(1, 1) >= 180 && Omega (1, 1) < 360 % Determining which value in the matrix is the larger one
        angle = Omega(1, 1);
    else 
        angle = Omega(1, 2);
    end
end
fprintf('%c = %.2f%c  //  %s --> %d%c \x2264 %c < %d%c\n', char(937), angle, char(176), type, angle_1, char(176), char(937), angle_2, char(176)); % Printing the value for RAAN and the reasoning 

% Argument of perigee
if e(1, 3) > 0 % Determining AoP based on value of the z-comp of eccentricity 
    type = 'e_z > 0';
    angle_1 = 0;
    angle_2 = 180;
    if omega(1, 1) >= 0 && omega (1, 1) < 180  % Determining which value in the matrix is the smaller one
        angle = omega(1, 1);
    else 
        angle = omega(1, 2);
    end 
else 
    type = 'e_z < 0';
    angle_1 = 180;
    angle_2 = 360;
    if omega(1, 1) >= 180 && omega (1, 1) < 360 % Determining which value in the matrix is the larger one
        angle = omega(1, 1);
    else 
        angle = omega(1, 2);
    end
end
fprintf('%c = %.2f%c  //  %s --> %d%c \x2264 %c < %d%c\n', char(969), angle, char(176), type, angle_1, char(176), char(969), angle_2, char(176)); % Printing the value for AoP and the reasoning

% True Anomaly
if v_r >= 0 % Determining true anomaly based on the value of v_r
    type = 'V_r > 0';
    effect = 'flying away from pergiee';
    angle_1 = 0;
    angle_2 = 180;
    if theta(1, 1) >= 0 && theta (1, 1) < 180 % Determining which value in the matrix is the smaller one
        angle = theta(1, 1);
    else 
        angle = theta(1, 2);
    end 
else 
    type = 'V_r < 0';
    effect = 'flying towards perigee';
    angle_1 = 180;
    angle_2 = 360;
    if theta(1, 1) >= 180 && theta (1, 1) < 360 % Determining which value in the matrix is the larger one
        angle = theta(1, 1);
    else 
        angle = theta(1, 2);
    end
end
fprintf('%c = %.2f%c  //  %s --> %d%c \x2264 %c < %d%c  //  Satellite is %s\n\n', char(952), angle, char(176), type, angle_1, char(176), char(952), angle_2, char(176), effect); %Printing the value of true anomaly and the reasoning

% Plotting Sketch
center = [0, 0, 0];
Num = 100;
t = linspace(0, 2*pi, Num);
x = center(1) + radius_E*cos(t);
y = center(2) + radius_E*sin(t);
z = center(3)*ones(size(x));
plot3(x, y, z, 'Color', 'blue'); % Plotting blue circle with radius = radius_e on XY plane
hold on;
x = center(3)*ones(size(x));
z = center(1) + radius_E*cos(t);
plot3(x, y, z, 'Color', 'blue'); % Plotting blue circle with radius = radius_e on YZ plane
hold on;
x = center(2) + radius_E*sin(t);
y = center(3)*ones(size(x));
plot3(x, y, z, 'Color', 'blue'); % Plotting blue circle with radius = radius_e on XZ plane
hold on;

% Finding greatest size vector
r_mags = [R_mag, V_mag]; % if you want to add more parameters, unit vectors will scale accordingly
r_max = max(r_mags);

% Define and scaling the unit vectors
i = [1, 0, 0]*r_max*1.1;
j = [0, 1, 0]*r_max*1.1;
k = [0, 0, 1]*r_max*1.1;

% Plotting the Orbital parameters in Geocentric Reference frame
scale = 4000;
% Position
quiver3(0, 0, 0, R(1), R(2), R(3), 'r', 'LineWidth', 2); % Plotting position vector
text(R(1), R(2), R(3), 'R', 'FontSize', 10); % Adding position title

% Velocity
quiver3(R(1), R(2), R(3), V(1)*scale, V(2)*scale/10, V(3)*scale, 'Color', 'green', 'LineWidth', 2); % Plotting velocity vector but scaled up for visual effect
text(R(1)+V(1)*scale, R(2)+V(2)*scale/10, R(3)+V(3)*scale,  'V', 'FontSize', 10); % Adding velocity title

% Displaying unit vectors
quiver3(0, 0, 0, i(1), i(2), i(3), 'k', 'LineWidth', 2); % Plotting unit vector I
text(i(1), i(2), i(3), 'I', 'FontSize', 12); % Adding unit vector I title
quiver3(0, 0, 0, j(1), j(2), j(3), 'k', 'LineWidth', 2); % Plotting unit vector J
text(j(1), j(2), j(3), 'J', 'FontSize', 12); % Adding unit vector J title
quiver3(0, 0, 0, k(1), k(2), k(3), 'k', 'LineWidth', 2); % Plotting unit vector K
text(k(1), k(2), k(3), 'K', 'FontSize', 12); % Adding unit vector K title

text(center(1), center(2), center(3), 'Earth', 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); % Adding Earth title
axis equal;
grid on;
xlabel('X'); % x-axis label
ylabel('Y'); % y-axis label
zlabel('Z'); % z-axis label
title('Position with respect the Geocentric reference frame') % title for entire plot

%% Part B - deltaV_total
r_p1 = (h_mag^2/mu_earth)*(1/(1 + e_mag)); % Perigee position for departure [km]

V_A1 = h_mag/r_p1; % V1 for departure from perigee position of ellipical parking orbit [km/s]
V_A2 = sqrt(mu_earth/r_p1); % Satellite circular parking orbit velocity [km/s]
deltaV_1 = V_A2 - V_A1; % Delta V required to go from elliptical orbit to circular parking orbit [km/s]

V_B1 = V_A2; % Satellite circular parking orbit velocity [km/s]
v_infi = sqrt(mu_sun/R_E)*(sqrt((2*R_J)/(R_E + R_J))-1); % Hyperbolic excess speed at Earth [km/s]
V_B2 = sqrt(v_infi^2 + ((2*mu_earth)/r_p1)); % Satellite helocentric velocity at perigee of Hohmann transfer [km/s]
deltaV_2 = V_B2 - V_B1; % Delta V required to go from Earth circular parking orbit to perigee of Hohmann transfer orbit [km/s]

v_infi = sqrt(mu_sun/R_J)*(sqrt((2*R_E)/(R_E + R_J))-1); % Hyperbolic excess speed at Jupiter [km/s]
V_C1 = sqrt(v_infi^2 + ((2*mu_jupiter)/(radius_J+200))); % Satellite helocentric velocity at apogee of Hohmann transfer [km/s]
V_C2 = sqrt(mu_jupiter/(radius_J+200)); % Satellite circular parking orbit velocity [km/s]
deltaV_3 = V_C2 - V_C1; % Delta V required to go from apogee of Hohmann transfer orbit to Jupiter parking orbit [km/s]

deltaV_total = abs(deltaV_3) + abs(deltaV_2) + abs(deltaV_1); % Total delta V required for one way mission [km/s]
fprintf('Part B: Total Delta V\n') % Printing Part title
deltaVs = [deltaV_1, deltaV_2, deltaV_3]; % All three delta V's in matrix form
for m = 1:3 % Printing all three delta V's
        if deltaVs(1, m) < 0
            type = '(Slowing down)';
        else
            type = '(Speeding up)';
        end
        fprintf('%cV_%d = %.2f km/s %s\n', char(916), m, deltaVs(1, m), type);
end
fprintf('%cV_Total = %.2f km/s\n\n', char(916), deltaV_total); % Printing total delta V

%% Part C - Hohmann Transfer orbit
a2 = (R_E + R_J)/2; % Semi-major axis of Hohmann transfer orbit [km]
e_hohmann = (R_J - R_E)/(R_J + R_E); % Eccentricity of Hohmann Transfer orbit
t_hohmann = ((2*pi()/sqrt(mu_sun))*(a2^(3/2)))/31536000; % Period of Hohmann transfer oribit [Earth years]

fprintf('Part C: Hohmann Transfer parameters\n')
fprintf('a_Hohmann = %.2e km\n', a2);
fprintf('e = %.5f\n', e_hohmann);
fprintf('t_Hohmann = %.2f Earth years\n\n', t_hohmann);

%% Part D - Phase angles
n1 = (2*pi())/t_E; % Mean motion of Earth's orbit [rad/year]
n2 = (2*pi())/t_J; % Mean motion of Jupiter's orbit [rad/year]
t_12 = t_hohmann/2; % Time it takes to traverse half of the Hohmann transfer ellipse [sec]
phi_i = pi() - n2*t_12; % Inital phase angle [rad]
phi_f = pi() - n1*t_12; % Final phase angle [rad]
T_syn = (t_E*t_J)/abs(t_E-t_J); % Synodic Period [sec]
if n1 > n2 % Determining if n1 > n2
    N = 0;
    t_wait = (-2*phi_f - 2*pi()*N)/(n2 - n1);
    while t_wait < 0
        N = N + 1; % Incrementing N
        t_wait = (-2*phi_f - 2*pi()*N)/(n2 - n1); % Wait time [sec]
    end
    
else % Determining if n2 > n1
    N = 0;
    t_wait = (-2*phi_f + 2*pi()*N)/(n2 - n1);
    
    while t_wait < 0
        N = N + 1; % Incrementing N
        t_wait = (-2*phi_f + 2*pi()*N)/(n2 - n1); % Wait time [sec]
    end
end
t_totaltrip = 2*t_12 + t_wait; % Total mission time [sec]

fprintf('Part D: Phase Angles and Time\n'); % Printing Part title
fprintf('%c_0 = %.3f rad \n', char(934), phi_i); % Printing intial phase angle
fprintf('%c_f = %.3f rad \n', char(934), phi_f); % Printing final phase angle
fprintf('T_synodic = %.2f Earth years \n', T_syn); % Printing synodic period
fprintf('t_12 (time of flight) = %.2f Earth years\n', t_12); % Printing time of flight
fprintf('t_wait (time spent waiting in Jupiter parking orbit) = %.2f Earth years \n', t_wait); % Printing wait time
fprintf('Total mission time = %.2f Earth years\n\n', t_totaltrip); % Printing total mission time

% Define the radius of the planets
r1 = 0.2;
r2 = 0.1;
r4 = 2;
r5 = 5;
lim = 5.5;

% Create a new figure
figure('Position', [100 100 1200 400]);
theta = linspace(0, 2*pi, 100);
subplot(1, 2, 1);

% Plot the EARTH    
center = [r4, 0];
x = r2*cos(theta) + center(1);
y = r2*sin(theta) + center(2);
plot(x, y, 'b');
hold on;

% Plotting EARTH - ORBIT
center = [0, 0];
x = r4*cos(theta) + center(1);
y = r4*sin(theta) + center(2);
plot(x, y, 'b');
hold on;

line([0, r5], [0, 0], 'Color', 'k')
hold on;

% Plot the parking orbit
% center = [2, 0];
% x = a*cos(theta) + center(1);
% y = b*sin(theta) + center(2);
% plot(x, y, 'k');
% hold on;

% Plotting JUPITER
center = [r5*cos(phi_i), r5*sin(phi_i)];
x = r2*cos(theta) + center(1);
y = r2*sin(theta) + center(2);
plot(x, y, 'r');
hold on;

% Plotting JUPITER - ORBIT
center = [0, 0];
x = r5*cos(theta) + center(1);
y = r5*sin(theta) + center(2);
plot(x, y, 'r');

line([0, r5*cos(phi_i)], [0, r5*sin(phi_i)], 'Color', 'k')
hold on;

% Adding arc to indicate phi_i  
r = (r4 + r5)/2;  % radius of the arc
arctheta = linspace(0, phi_i, 100);  % angle range for the arc
x_arc = r*cos(arctheta);  % x-coordinates of the arc
y_arc = r*sin(arctheta);  % y-coordinates of the arc
plot(x_arc, y_arc, 'k')  % plot the arc 
text((r*cos(phi_i/2))-0.1, (r*sin(phi_i/2))-0.1, 'Φ_0', 'BackgroundColor', 'w', 'EdgeColor', 'w', 'Margin', 1)
hold on;

% Plot the SUN
center = [0, 0];
x = r1*cos(theta) + center(1);
y = r1*sin(theta) + center(2);
plot(x, y, 'Color', [1, 0.5, 0]);
hold on;

% Set the axis limits to show both circles
xlim([-lim, lim]);
ylim([-lim, lim]);
axis equal;
% Add a title to the plot
title('Planet positions upon satellite departure (not to scale)');

subplot(1, 2, 2);
% Plotting EARTH
center = [r4*cos(t_12*(2*pi()/t_E)), r4*sin(t_12*(2*pi()/t_E))];
x = r2*cos(theta) + center(1);
y = r2*sin(theta) + center(2);
plot(x, y, 'b');
hold on;

% Plotting EARTH - ORBIT
center = [0, 0];
x = r4*cos(theta) + center(1);
y = r4*sin(theta) + center(2);
plot(x, y, 'b');
hold on;

line([0, r5*cos(t_12*(2*pi()/t_E))], [0, r5*sin(t_12*(2*pi()/t_E))], 'Color', 'k')
hold on;

% Plotting JUPITER
center = [r5*cos((t_12*(2*pi()/t_J)) + phi_i), r5*sin((t_12*(2*pi()/t_J)) + phi_i)];
x = r2*cos(theta) + center(1);
y = r2*sin(theta) + center(2);
plot(x, y, 'r');
hold on;

% Plotting JUPITER - ORBIT
center = [0, 0];
x = r5*cos(theta) + center(1);
y = r5*sin(theta) + center(2);
plot(x, y, 'r');

line([0, r5*cos((t_12*(2*pi()/t_J)) + phi_i)], [0, r5*sin((t_12*(2*pi()/t_J)) + phi_i)], 'Color', 'k')
hold on;

angle2_E = t_12*(2*pi()/t_E);
angle2_J = (t_12*(2*pi()/t_J)) + phi_i;
while angle2_E > 2*pi()
    angle2_E = angle2_E - 2*pi();
end
while angle2_J > 2*pi()
    angle2_J = angle2_J - 2*pi();
end
% Adding arc to indicate phi_f  
r = (r4 + r5)/2;  % radius of the arc
arctheta = linspace(angle2_E, angle2_J, 100);  % angle range for the arc
x_arc = r*cos(arctheta);  % x-coordinates of the arc
y_arc = r*sin(arctheta);  % y-coordinates of the arc
plot(x_arc, y_arc, 'k')  % plot the arc 
if angle2_E > angle2_J
    text((r*cos(abs(angle2_E - angle2_J)/2 + angle2_J))-0.1, (r*sin(abs(angle2_E - angle2_J)/2 + angle2_J))-0.1, 'Φ_f', 'BackgroundColor', 'w', 'EdgeColor', 'w', 'Margin', 1)
else
    text((r*cos(abs(angle2_E - angle2_J)/2 + angle2_E))-0.1, (r*sin(abs(angle2_E - angle2_J)/2 + angle2_E))-0.1, 'Φ_f', 'BackgroundColor', 'w', 'EdgeColor', 'w', 'Margin', 1)
end
    hold on;

% Plot the SUN
center = [0, 0];
x = r1*cos(theta) + center(1);
y = r1*sin(theta) + center(2);
plot(x, y, 'Color', [1, 0.5, 0]);
hold on;

xlim([-lim, lim]);
ylim([-lim, lim]);
axis equal;
title("Planet positions upon satellite arrival (not to scale)")
%% Part E - Propellant Choice
fprintf('Part E: Propellant Choice\n')
m_ratio_solid = (1-exp(-deltaV_2/(I_sp_solid*g_0)))*100;
m_ratio_liquidH = (1-exp(-deltaV_2/(I_sp_liquidH*g_0)))*100;
if m_ratio_solid > m_ratio_liquidH
    fuel = 'Liquid hydrogen';
    lesser = 'LH';
    greater = 'S';
else 
    fuel = 'Solid';
    lesser = 'S';
    greater = 'LH';
end
fprintf('%cm/m0_S = %.2f%% (I_sp = %d sec)\n', char(916), m_ratio_solid, I_sp_solid);
fprintf('%cm/m0_LH = %.2f%% (I_sp = %d sec)\n', char(916), m_ratio_liquidH, I_sp_liquidH);
fprintf('Optimal propellant type for mission: %s because %cm/m0_%s < %cm/m0_%s\n', fuel, char(916), lesser, char(916), greater)
% Can't finish bc I have no idea what I'm doing with this