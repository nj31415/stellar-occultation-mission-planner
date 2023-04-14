%% Nicholas Jones - njones31@vt.edu
% Satellite
% The satellite class contains the data pertinent to a satellite as defined
% in STK, as well as the data needed for occultation analysis.

classdef Satellite < handle
    properties (Constant)
        mu = 398600.4418; % km^3 s^-2 - Earth gravitational parameter
    end
    
    properties (SetAccess = public)
        name

        % Initial Classical Orbital Elements
        a_init
        e_init
        i_init
        lan_init
        omega_init
        nu_init

        % Spacecraft aerodynamic properties
        drag_coeff
        area_mass_ratio

        % Ephemeris data. Should be column vectors
        time_stamp = [];
        x_icrf = [];
        y_icrf = [];
        z_icrf = [];
        vx_icrf = [];
        vy_icrf = [];
        vz_icrf = [];
    end

    properties (SetAccess = private)
        % Occultation data. Setting operations provided by add and remove
        % methods.
        occ_data = Occ_Event.empty;
    end

    properties (Access = private)
        dt_helper   % Datetime_helper object
    end
    
    methods
        %% Constructor for a Satellite object. All values are referenced to
        % STK's TrueOfDate coordinate system. Throws errors when
        % encountering data that violates conditions noted in Inputs below.
        %
        % Inputs:
        % scenario_start    : String, formatted string representing the
        %                     scenario start time, follows format specified
        %                     in Datetime_Helper constructor.
        % name              : String, identifier for the satellite object.
        %                     Must follow STK naming conventions.
        % a_init            : Float, initial semimajor axis in km.
        % e_init            : Float, initial eccentricity.
        % i_init            : Float, initial inclination in degrees. Must
        %                     be between 0 and 180 degrees.
        % lan_init          : Float, initial Longitude of Ascending Node
        %                     in degrees. Must be between 0 and 360 degrees
        % omega_init        : Float, initial argument of perigee in
        %                     degrees. Must be between 0 and 360 degrees.
        % nu_init           : Float, initial true anomaly in degrees. Must
        %                     be between 0 and 360 degrees.
        % drag_coeff        : Float, satellite drag coefficient for
        %                     lifetime analysis.
        % area_mass_ratio   : Float, satellite area to mass ratio for
        %                     lifetime analysis
        function obj = Satellite(scenario_start, name, a_init, e_init, ...
            i_init, lan_init, omega_init, nu_init, drag_coeff, ...
            area_mass_ratio)
            obj.dt_helper = Datetime_Helper(scenario_start);
            obj.name = name;
            obj.a_init = a_init;
            obj.e_init = e_init;
            obj.i_init = i_init;
            obj.lan_init = lan_init;
            obj.omega_init = omega_init;
            obj.nu_init = nu_init;
            obj.drag_coeff = drag_coeff;
            obj.area_mass_ratio = area_mass_ratio;
        end
        
        %% Get the satellite name.
        % Outputs:
        % name  : String, the satellite name
        function name = get.name(obj)
            name = obj.name;
        end

        %% Set the satellite name.
        % inputs:
        % name  : String, identifier of the satellite object. Must follow
        %         STK naming conventions.
        function set.name(obj, name)
            if ~contains(name, ' ')
                obj.name = name;
            else
                error('Satellite names cannot have spaces');
            end
        end

        %% Get the satellite initial semi-major axis in km
        % Outputs:
        % a_init    : Float, initial semimajor axis in km
        function a_init = get.a_init(obj)
            a_init = obj.a_init;
        end

        %% Set the satellite initial semi-major axis in km
        % Inputs:
        % a_init    : Float, initial semimajor axis in km
        function set.a_init(obj, a_init)
            obj.a_init = a_init;
        end

        %% Get the satellite initial eccentricity
        % Outputs:
        % e_init    : Float, initial eccentricity
        function e_init = get.e_init(obj)
            e_init = obj.e_init;
        end

        %% Set the satellite initial eccentricity
        % Inputs:
        % e_init    : Float, initial eccentricity
        function set.e_init(obj, e_init)
            obj.e_init = e_init;
        end

        %% Get the satellite initial inclination in degrees.
        % Outputs:
        % i_init    : Float, initial inclination in degrees
        function i_init = get.i_init(obj)
            i_init = obj.i_init;
        end

        %% Set the satellite initial inclination in degrees between 0 and
        % 180 degrees.
        % Inputs:
        % i_init    : Float, initial inclination in degrees, between 0 and
        %             180 degrees
        function set.i_init(obj, i_init)
            if i_init >= 0 && i_init <= 180
                obj.i_init = i_init;
            else
                error(['Satellite::set.i_init: Inclination falls'...
                    'outisde the acceptable bounds']);
            end
        end

        %% Get the satellite initial Longitude of Ascending Node in
        % degrees.
        % Outputs:
        % lan_init  : Float, Longitude of Ascending Node in degrees.
        function lan_init = get.lan_init(obj)
            lan_init = obj.lan_init;
        end

        %% Set the satellite initial Longitude of Ascending Node in
        % degrees. Must be between 0 and 360 degrees.
        % Inputs:
        % lan_init  : Float, Longitude of Ascending Node in degrees,
        %             between 0 and 360 degrees.
        function set.lan_init(obj, lan_init)
            if lan_init >= 0 && lan_init <= 360
                obj.lan_init = lan_init;
            else
                error(['Satellite::set.lan_init: Longitude of Ascending'...
                    'Node falls outisde the acceptable bounds']);
            end
        end

        %% Get the satellite initial argument of periapsis in degrees.
        % Outputs:
        % omega_init    : Float, initial argument of periapsis in degrees.
        function omega_init = get.omega_init(obj)
            omega_init = obj.omega_init;
        end

        %% Set the satellite initial argument of periapsis in degrees. Must
        % be between 0 and 360 degrees.
        % Inputs:
        % omega_init    : Float, initial argument of periapsis in degrees,
        %                 between 0 and 360 degrees.
        function set.omega_init(obj, omega_init)
            if omega_init >= 0 && omega_init <= 360
                obj.omega_init = omega_init;
            else
                error(['Satellite::set.omega_init: Argument of '...
                    'Periapsis falls outisde the acceptable bounds']);
            end
        end

        %% Get the satellite initial true anomaly in degrees.
        % Outputs:
        % nu_init   : Float, initial true anomaly in degrees.
        function nu_init = get.nu_init(obj)
            nu_init = obj.nu_init;
        end

        %% Set the satellite initial true anomaly in degrees. Must be
        % between 0 and 360 degrees.
        % Inputs:
        % nu_init   : Float, initial true anomaly, between 0 and 360
        %             degrees.
        function set.nu_init(obj, nu_init)
            if nu_init >= 0 && nu_init <= 360
                obj.nu_init = nu_init;
            else
                error(['Satellite::set.nu_init: True Anomaly falls'...
                    'outisde the acceptable bounds']);
            end
        end

        %% Get the satellite drag coefficient.
        % Outputs:
        % drag_coeff    : Float, Drag coefficient of the satellite
        function drag_coeff = get.drag_coeff(obj)
            drag_coeff = obj.drag_coeff;
        end

        %% Set the satellite drag coefficient
        % Inputs:
        % drag_coeff    : Float, drag coefficient of the satellite
        function set.drag_coeff(obj, drag_coeff)
            obj.drag_coeff = drag_coeff;
        end

        %% Get the satellite drag coefficient.
        % Outputs:
        % drag_coeff    : Float, Drag coefficient of the satellite
        function area_mass_ratio = get.area_mass_ratio(obj)
            area_mass_ratio = obj.area_mass_ratio;
        end

        %% Set the satellite drag coefficient
        % Inputs:
        % drag_coeff    : Float, drag coefficient of the satellite
        function set.area_mass_ratio(obj, area_mass_ratio)
            obj.area_mass_ratio = area_mass_ratio;
        end

        %% Get the satellite time stamps for ephemeris data in EpSec.
        % Outputs:
        % time_stamp    : Float vector, time stamps in EpSec as a column
        %                 vector.
        function time_stamp = get.time_stamp(obj)
            time_stamp = obj.time_stamp;
        end

        %% Get the last propagated time stamp in EpSec
        % Outputs:
        % time_stamp_end    : Float, final time stamp in EpSec
        function time_stamp_end = lifetime(obj)
            time_stamp_end = obj.time_stamp(end);
        end

        %% Set the satellite time stamps for ephemeris data in EpSec. Must
        % be the same length as the other (position, velocity) ephemeris data.
        % Inputs:
        % time_stamp  : Float vector, time stamps in EpSec as a column
        %               vector.
        function set.time_stamp(obj, time_stamp)
            obj.time_stamp = reshape(time_stamp, [], 1);
        end
        
        %% Get the satellite ICRF X position ephemeris data in km.
        % Outputs:
        % x_icrf    : Float vector, X ICRF positions in km as a column
        %             vector.
        function x_icrf = get.x_icrf(obj)
            x_icrf = obj.x_icrf;
        end

        %% Set the satellite ICRF X position ephemeris data in km. Must be
        % the same length as the other (position, velocity) ephemeris data.
        % Inputs:
        % x_icrf    : Float vector, X ICRF positions in km as a column
        %             vector.
        function set.x_icrf(obj, x_icrf)
            obj.x_icrf = reshape(x_icrf, [], 1);
        end

        %% Get the satellite ICRF Y position ephemeris data in km.
        % Outputs:
        % y_icrf    : Float vector, Y ICRF positions in km as a column
        %             vector.
        function y_icrf = get.y_icrf(obj)
            y_icrf = obj.y_icrf;
        end

        %% Set the satellite ICRF Y position ephemeris data in km. Must be
        % the same length as the other (position, velocity) ephemeris data.
        % Inputs:
        % y_icrf    : Float vector, Y ICRF positions in km as a column
        %             vector.
        function set.y_icrf(obj, y_icrf)
            obj.y_icrf = reshape(y_icrf, [], 1);
        end

        %% Get the satellite ICRF Z position ephemeris data in km.
        % Outputs:
        % z_icrf    : Float vector, Z ICRF positions in km as a column
        %             vector.
        function z_icrf = get.z_icrf(obj)
            z_icrf = obj.z_icrf;
        end

        %% Set the satellite ICRF Z position ephemeris data in km. Must be
        % the same length as the other (position, velocity) ephemeris data.
        % Inputs:
        % z_icrf    : Float vector, Z ICRF positions in km as a column
        %             vector.
        function set.z_icrf(obj, z_icrf)
            obj.z_icrf = reshape(z_icrf, [], 1);
        end

        %% Get the satellite ICRF positions ephemeris data in km using a
        % single matrix.
        % Inputs:
        % ind           : Index vector, the desired indices to retrieve. If
        %                 unspecified, retrieves all data.
        % Outputs:
        % icrf  : Float matrix, 1st column is X, 2nd is Y, and 3rd is Z
        %         ICRF position ephemeris data in km.
        function icrf = get_icrf_position(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                icrf = [obj.x_icrf(ind), obj.y_icrf(ind), obj.z_icrf(ind)];
            else
                icrf = [obj.x_icrf, obj.y_icrf, obj.z_icrf];
            end
        end

        %% Set the satellite ICRF position ephemeris data in km using a 
        % single function. Inputs should be column vectors.
        % Inputs:
        % x_icrf    : Float vector, X ICRF positions in km as a column
        %             vector
        % y_icrf    : Float vector, Y ICRF positions in km as a column
        %             vector
        % z_icrf    : Float vector, Z ICRF positions in km as a column
        %             vector
        function set_icrf_position(obj, x_icrf, y_icrf, z_icrf)
            if length(x_icrf) == length(y_icrf) && ...
                    length(x_icrf) == length(z_icrf)
                obj.x_icrf = x_icrf;
                obj.y_icrf = y_icrf;
                obj.z_icrf = z_icrf;
            else
                error(['Satellite::set_icrf_position: Size of input '...
                    'vectors are mismatched']);
            end
        end

        %% Get the satellite ICRF X velocity ephemeris data in km s^-1.
        % Outputs:
        % vx_icrf    : Float vector, X ICRF velocities in km s^-1 as a
        %              column vector.
        function vx_icrf = get.vx_icrf(obj)
            vx_icrf = obj.vx_icrf;
        end

        %% Set the satellite ICRF X velocity ephemeris data in km. Must be
        % the same length as the other (position, velocity) ephemeris data.
        % Inputs:
        % vx_icrf    : Float vector, X ICRF velocities in km s^-1.
        function set.vx_icrf(obj, vx_icrf)
            obj.vx_icrf = reshape(vx_icrf, [], 1);
        end

        %% Get the satellite ICRF Y velocity ephemeris data in km s^-1.
        % Outputs:
        % vy_icrf    : Float vector, Y ICRF velocities in km s^-1 as a
        %              column vector.
        function vy_icrf = get.vy_icrf(obj)
            vy_icrf = obj.vy_icrf;
        end

        %% Set the satellite ICRF Y velocity ephemeris data in km. Must be
        % the same length as the other (position, velocity) ephemeris data.
        % Inputs:
        % vy_icrf    : Float vector, Y ICRF velocities in km s^-1 as a
        %              column vector.
        function set.vy_icrf(obj, vy_icrf)
            obj.vy_icrf = reshape(vy_icrf, [], 1);
        end

        %% Get the satellite ICRF Z velocity ephemeris data in km s^-1.
        % Outputs:
        % vz_icrf    : Float vector, Z ICRF velocities in km s^-1 as a
        %              column vector.
        function vz_icrf = get.vz_icrf(obj)
            vz_icrf = obj.vz_icrf;
        end

        %% Set the satellite ICRF Z velocity ephemeris data in km. Must be
        % the same length as the other (position, velocity) ephemeris data.
        % Inputs:
        % vz_icrf    : Float vector, Z ICRF velocities in km s^-1 as a
        %              column vector.
        function set.vz_icrf(obj, vz_icrf)
            obj.vz_icrf = reshape(vz_icrf, [], 1);
        end

        %% Get the satellite ICRF velocities ephemeris data in km s^-1
        % using a single matrix.
        % Inputs:
        % ind           : Index vector, the desired indices to retrieve. If
        %                 unspecified, retrieves all data.
        % Outputs:
        % v_icrf    : Float matrix, 1st column is X, 2nd is Y, and 3rd is Z
        %             ICRF velocity ephemeris data in km s^-1.
        function v_icrf = get_icrf_velocity(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                v_icrf = [obj.vx_icrf(ind), obj.vy_icrf(ind), ...
                    obj.vz_icrf(ind)];
            else
                v_icrf = [obj.vx_icrf(:), obj.vy_icrf(:), ...
                    obj.vz_icrf(:)];
            end
        end

        %% Set the satellite ICRF velocities ephemeris data in km s^-1
        % using a single function. Inputs should be column vectors.
        % Inputs:
        % vx_icrf   : Float vector, X ICRF velocities in km s^-1 as a
        %             column vector
        % vy_icrf   : Float vector, Y ICRF velocities in km s^-1 as a
        %             column vector
        % vz_icrf   : Float vector, Z ICRF velocities in km s^-1 as a
        %             column vector
        function set_icrf_velocity(obj, vx_icrf, vy_icrf, vz_icrf)
            if length(vx_icrf) == length(vy_icrf) &&...
                    length(vx_icrf) == length(vz_icrf)
                obj.vx_icrf = vx_icrf;
                obj.vy_icrf = vy_icrf;
                obj.vz_icrf = vz_icrf;
            else
                error(['Satellite::set_icrf_velocity: Size of input '...
                    'vectors are mismatched']);
            end
        end

        %% Get the occultatation data vector for a satellite. An entry is
        % composed of an Occ_Event object.
        % Outputs:
        % occ_data  : Occ_Event vector, occultation data for the satellite.
        function occ_data = get.occ_data(obj)
            occ_data = obj.occ_data;
        end

        %% Get the occultation data for a satellite for a specified star.
        % Inputs:
        % star_name : String array, name(s) of the star to retrieve data
        %             for.
        function [occ_event] = search_occ_data_star(obj, star_name)
            search_idx = matches({obj.occ_data.star_name}, star_name);
            occ_event = obj.occ_data(search_idx);
        end

        %% Add occultation data to the occultation data vector for a
        % satellite.
        % Inputs:
        % occ_event : Occ_Event, an Occ_Event object associated with a
        %             star.
        function add_occ_data(obj, occ_event)
            for i = 1 : length(occ_event)
                    obj.occ_data(1, end + 1) = occ_event(i);
            end
        end

        %% Remove occultation data from the occultation data vector for a
        % satellite for a specified star.
        % Inputs:
        % star_name : String array, name(s) of the stars to remove data
        %             from the occultation data vector.
        function remove_occ_data(obj, star_name)
            search_idx = matches({obj.occ_data.star_name}, star_name);
            obj.occ_data(search_idx) = [];
        end

        %% Clear the occultation data vector for the satellite.
        function clear_occ_data(obj)
            obj.occ_data = [];
        end

        %% Retrieve logical indices for satellite ephemeris data within
        % a given time interval. Searches are inclusive of the limits. Can
        % be used as an input for the getter functions within this class.
        % Inputs:
        % time_lim  : Float vector, specification of the minimum and
        %             maximum times in EpSec. [min_time max_time]
        % Outputs:
        % ind_time  : Logical index array, true when the tangent point
        %             occurs within the time period of interest.
        function ind_time = search_time(obj, time_lim)
            ind_time = obj.time_stamp(:, 1) >= time_lim(1) & ...
                obj.time_stamp(:, 1) <= time_lim(2);
        end

        %% Calculate the classical orbital elements from the satellite
        % ephemeris data
        % Inputs:
        % ind           : Index vector, the desired indices to retrieve.
        % Outputs:
        % a     : Float vector, semimajor axis in km
        % e     : Float vector, eccentricity in degrees
        % i     : Float vector, inclination in degrees
        % RAAN  : Float vector, right ascension of ascending node in
        %         degrees
        % omega : Float vector, argument of periapsis in degrees
        % nu    : Float vector, true anamoly in degrees
        function [a, e, i, RAAN, omega, nu] = ...
                ephem_2_classical_elem(obj, ind)
            r = obj.get_icrf_position(ind);
            v = obj.get_icrf_velocity(ind);

            h = cross(r, v, 2);
            n = cross([zeros(size(h, 1), 2) ones(size(h, 1), 1)], h, 2);
            
            v_norm = vecnorm(v, 2, 2);
            r_norm = vecnorm(r, 2, 2);

            e_vec = ((v_norm.^2 - (obj.mu ./ r_norm)) .* r - ...
                dot(r, v, 2) .* v) ./ obj.mu;

            e = vecnorm(e_vec, 2, 2);

            epsilon = (v_norm.^2 / 2) - obj.mu ./ r_norm;

            a = -obj.mu ./ (2 * epsilon);

            i = acosd(h(:, 3) ./ vecnorm(h, 2, 2));

            RAAN = acosd(n(:, 1) ./ vecnorm(n, 2, 2));
            RAAN(n(:, 2) < 0) = 360 - RAAN(n(:, 2) < 0);

            omega = acosd(dot(n, e_vec, 2) ./ (vecnorm(n, 2, 2) .* ...
                vecnorm(e_vec, 2, 2)));
            omega(e_vec(:, 3) < 0) = 360 - omega(e_vec(:, 3) < 0);

            nu = acosd(dot(e_vec, r, 2) ./ (vecnorm(e_vec, 2, 2) .* ...
                vecnorm(r, 2, 2)));
            nu(dot(r, v, 2) < 0) = 360 - nu(dot(r, v, 2) < 0);
        end
    end
end