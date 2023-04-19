%% Nicholas Jones - njones31@vt.edu
% Plot change in inclination with semimajor axis for Sun Synchronous
% Orbits.
clear;
clc;

R_E = physconst('EarthRadius') * (10^-3); % km
a = R_E + linspace(0, 1000, 100); % km - semimajor axis

e = 0; % Eccentricity;

RAAN_dot_SS = (2 * pi / 365.2421897) * (1 / 24) * (1 / 3600); % rad sec^-1
J2 = 0.0010826267;  % J2 Constant
mu = 398600.4418; % km^3 sec^-2 - Earth Gravitational Parameter

inc = @(a) acosd((-2 * a.^(7 / 2) * RAAN_dot_SS * (1 - (e^2))^2) ...
    / (3 * R_E^2 * J2 * sqrt(mu)));

figure();
plot(a - R_E, inc(a), 'k*');
xlabel('Initial Altitude (km)');
ylabel('Inclination (deg)');
title('Satellite Altitude vs Inclination for Sun-Synchronous Orbits');
xline(350);
xline(400);
xline(450);