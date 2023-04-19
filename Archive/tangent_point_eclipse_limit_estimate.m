%% Nicholas Jones - njones31@vt.edu
% Estimate of tangent point eclipse limits
% Estimate the altitude limit at which a tangent point is still in eclipse,
% based on a spherical Earth with cylindrical shadow.
close all;
clear;
clc;

R_E = linspace(6357, 6378, 5)';
theta = linspace(60, 90, 31);

h = R_E .* (sqrt(1 + cotd(theta)) - 1);

plot(90 - theta, h);
yline(150);
leg = legend(num2str(R_E));
title(leg, 'Earth Radius (km)');
ylabel('Altitude Eclipse Limit (km)');
xlabel('Angle from Terminator (deg)');
title('Tangent Point Eclipse Limit for Different Values of Earth Radius');