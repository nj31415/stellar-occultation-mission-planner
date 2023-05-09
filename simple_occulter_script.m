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
    [rp_1, ~] = occult(r_s_1, r_ts_1(i, :) .* ones(size(r_s_1)));
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
    [rp_2, ~] = occult(r_s_2, r_ts_2(i, :) .* ones(size(r_s_2)));
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

%% Study 3: Altitude rate of change for two orbits
r1 = 350 + R_E;
r2 = 450 + R_E;
h = (70 : 150) + R_E;
v_ratio = sqrt(r2 / r1) * sin(asin(h / r1) ./ asin(h / r2));

figure();
plot(h - R_E, v_ratio, 'k*');
xlabel('Tangent Point Altitude, km');
ylabel('Altitude Rate of Change Ratio');

%% Study 4: Snapshots of orbit plane angle study
r_alt_4 = 450 + R_E;
nu_4 = (0 : 359)';

r_s_4 = sat_gen(r_alt_4, nu_4);

theta_4 = (180 + ([0 30 60 89]))';
lambda_4 = 0;
r_ts_4 = sphere2cart(1, lambda_4 * ones(size(theta_4)), theta_4);

for i = 1 : length(theta_4)
    figure();
    [rp_4, rs_f_4] = occult(r_s_4, r_ts_4(i, :) .* ones(size(r_s_4)));
    rp_4_s = cart2sphere(rp_4(:, 1), rp_4(:, 2), rp_4(:, 3));
    rp_4(rp_4_s(:, 1) < R_E | isnan(rp_4_s(:, 1)), :) = NaN;
    rs_f_4(rp_4_s(:, 1) < R_E | isnan(rp_4_s(:, 1)), :) = NaN;
    
    plot3(rs_f_4(:, 1), rs_f_4(:, 2), rs_f_4(:, 3), 'r*');
    hold on;
    plot3(rp_4(:, 1), rp_4(:, 2), rp_4(:, 3), 'b*');
    plot3(r_s_4(:, 1), r_s_4(:, 2), r_s_4(:, 3), 'c');

    % Plot Earth - sphere on figure
    [X_s, Y_s, Z_s] = sphere;
    X_s = X_s * R_E;
    Y_s = Y_s * R_E;
    Z_s = Z_s * R_E;
    surf(X_s, Y_s, Z_s);
    colormap summer;

    quiver3(0, 0, 0, r_ts_4(i, 1), r_ts_4(i, 2), r_ts_4(i, 3), ...
        R_E * 1.25, 'LineWidth', 1.5);

    xlabel('X, km');
    ylabel('Y, km');
    zlabel('Z, km');
    axis equal;
    grid on;
    xlim(1.25 * [-R_E R_E]);
    ylim(1.25 * [-R_E R_E]);
    zlim(1.25 * [0 R_E]);
end

%% Study 5: Snapshots of declination study
r_alt_5 = 450 + R_E;
nu_5 = (0 : 359)';

r_s_5 = sat_gen(r_alt_5, nu_5);

theta_5 = 180;
lambda_5 = [0 30 60 90]';
r_ts_5 = sphere2cart(1, lambda_5 * ones(size(theta_5)), theta_5);

for i = 1 : length(lambda_5)
    figure();
    [rp_5, rs_f_5] = occult(r_s_5, r_ts_5(i, :) .* ones(size(r_s_5)));
    rp_5_s = cart2sphere(rp_5(:, 1), rp_5(:, 2), rp_5(:, 3));
    rp_5(rp_5_s(:, 1) < R_E | isnan(rp_5_s(:, 1)), :) = NaN;
    rs_f_5(rp_5_s(:, 1) < R_E | isnan(rp_5_s(:, 1)), :) = NaN;
    
    plot3(rs_f_5(:, 1), rs_f_5(:, 2), rs_f_5(:, 3), 'r*');
    hold on;
    plot3(rp_5(:, 1), rp_5(:, 2), rp_5(:, 3), 'b*');
    plot3(r_s_5(:, 1), r_s_5(:, 2), r_s_5(:, 3), 'c');

    % Plot Earth - sphere on figure
    [X_s, Y_s, Z_s] = sphere;
    X_s = X_s * R_E;
    Y_s = Y_s * R_E;
    Z_s = Z_s * R_E;
    surf(X_s, Y_s, Z_s);
    colormap summer;

    quiver3(0, 0, 0, r_ts_5(i, 1), r_ts_5(i, 2), r_ts_5(i, 3), ...
        R_E * 1.25, 'LineWidth', 1.5);

    xlabel('X, km');
    ylabel('Y, km');
    zlabel('Z, km');
    axis equal;
    grid on;
    xlim(1.25 * [-R_E R_E]);
    ylim(1.25 * [-R_E R_E]);
    zlim(1.25 * [-R_E R_E]);
    view(0, 0);
end

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
function [rp, rs_f] = occult(rs, rts)
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
rp(dot_vec >= 0, :) = NaN;

rs_f = rs;
rs_f(dot_vec >= 0, :) = NaN;

end