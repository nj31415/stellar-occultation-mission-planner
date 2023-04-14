% Nicholas Jones - njones31@vt.edu
% Script for exploring stellar occultation through a simplified model.
close all;
clear;
clc;

% Define constants
R_E = physconst('EarthRadius') * 10^-3;

%% Study 1: Impact of angle between the target star and orbit plane on
% latitude smearing of the tangent point during an observation
r_alt_1 = 450 + R_E;
nu_1 = (0 : 359)';

r_s_1 = sat_gen(r_alt_1, nu_1);

theta_1 = (180 + (-0 : 90))';
lambda_1 = 0;
r_ts_1 = sphere2cart(1, lambda_1 * ones(size(theta_1)), theta_1);

figure();

for i = 1 : length(theta_1)
    rp_1 = occult(r_s_1, r_ts_1(i, :) .* ones(size(r_s_1)));
    rp_1 = cart2sphere(rp_1(:, 1), rp_1(:, 2), rp_1(:, 3));
    rp_1(rp_1(:, 1) < R_E, :) = NaN;

    scatter((theta_1(i) - 180) * ones(size(rp_1(:, 1))), rp_1(:, 3), ...
        [], rp_1(:, 1) - R_E, 'filled');
    hold on;
end

xlabel('Angle from Orbit Plane, degrees');
ylabel('Latitude of Tangent Point, degrees');
ylim([-90 90]);
c = colorbar;
c.Label.String = 'Tangent Point Altitude, km';
colormap jet;

%% Study 2: Impact of star declination on the potential tangent point
% latitudes
r_alt_2 = 450 + R_E;
nu_2 = (0 : 359)';

r_s_2 = sat_gen(r_alt_2, nu_2);

lambda_2 = (-90 : 90)';
r_ts_2 = sphere2cart(1, lambda_2, 180 * ones(size(lambda_2)));

figure();

for i = 1 : length(lambda_2)
    rp_2 = occult(r_s_2, r_ts_2(i, :) .* ones(size(r_s_2)));
    rp_2 = cart2sphere(rp_2(:, 1), rp_2(:, 2), rp_2(:, 3));
    rp_2(rp_2(:, 1) < R_E, :) = NaN;

    plot(lambda_2(i) * ones(size(rp_2(:, 1))), rp_2(:, 3), 'k*');
    hold on;
end

xlabel('Star Declination, degrees');
ylabel('Latitude of Tangent Point, degrees');
ylim([-90 90]);
xlim([-90 90]);
yline(66.56361, '--k', 'Label', 'Arctic Circle', 'LineWidth', 1);
yline(-66.56361, '--k', 'Label', 'Antarctic Circle', 'LineWidth', 1);

% Function to generate satellite cartesian coordinates from an altitude and
% true anomaly.
function r_sat = sat_gen(r_alt, nu)
r_sat = [r_alt .* cosd(nu), zeros(size(nu)), r_alt .* sind(nu)];
end

% Function to change from cartesian to spherical coordinates
function r_cart = cart2sphere(x, y, z)
r_cart = [sqrt(x.^2 + y.^2 + z.^2), ...
    atan2d(y, x), ...
    atan2d(z, sqrt(x.^2 + y.^2))];
end

% Function to change from spherical to cartesian coordinates
function r_sphere = sphere2cart(alt, lambda, theta)
r_sphere = [alt .* cosd(lambda) .* cosd(theta), ...
    alt .* cosd(lambda) .* sind(theta), ...
    alt .* sind(lambda)];
end

% Function to calculate the location of the tangent point. Based on
% occult_1 function of Occulter.m
function rp = occult(rs, rts)
dot_vec = dot(rs, rts, 2);

cross_vec = cross(rs, rts, 2);
norm_cross_vec = vecnorm(cross_vec, 2, 2);

phi = atan2(norm_cross_vec, dot_vec);
theta = phi - (pi / 2);

point_alt = vecnorm(rs, 1, 2) .* cos(theta);

r_hat = cross_vec ./ norm_cross_vec;

v = rs ./ vecnorm(rs, 1, 2);

rp = (((1 - cos(theta)) .* dot(v, r_hat, 2)) .* r_hat + cos(theta) .* v ...
    + sin(theta) .* cross(r_hat, v, 2));
rp = point_alt .* rp;

% Filter invalid occultations
rp(dot_vec >= 0) = NaN;

end