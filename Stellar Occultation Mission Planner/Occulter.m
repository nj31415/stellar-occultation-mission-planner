%% Nicholas Jonse - njones31@vt.edu
% Occulter
% The Occulter class calculates the occultation geometry and populates an
% array of Occ_Events for a provided star and satellite, along
% with access data from STK Driver. Currently only a single function, but
% can be expanded with other methods of calculating the occultation
% geometry in the future with different inputs.

classdef Occulter < handle
    properties (Access = private)
        dt_helper       % Datetime_Helper object
        scenario_start  % start time placeholder for creating
                        % Datetime_Helper objects in the Occ_Event objects
    end

    methods
        %% Constructor for an Occulter object. Used to access methods for
        % calculating the occultation geometry.
        % Inputs:
        % scenario_start    : String, formatted string representing the
        %                     scenario start time, follows format specified
        %                     in Datetime_Helper constructor.
        function obj = Occulter(scenario_start)
            obj.dt_helper = Datetime_Helper(scenario_start);
            obj.scenario_start = scenario_start;
        end
        
        %% Function to create an Occ_Event object for a star and
        % corresponding satellite. Only considers the occultation geometry,
        % does not consider refraction.
        % Inputs:
        % star_name : String, name of the star being occulted
        % time      : Float vector, time stamps in EpSec during valid
        %             access periods between the star and the satellite.
        % r_sat     : Float vector, ICRF satellite position in km during
        %             valid access periods between the star and the
        %             satellite.
        % r_star    : Float vector, ICRF star position (normalized) during
        %             valid access periods between the star and the
        %             satellite.
        % time_step : Float, time step of the simulation in sec
        % Outputs:
        % occ_events: Occ_Event vector, completed Occ_Events array for the 
        %             satellite and star.
        function occ_events = occult_1(obj, star_name, time, r_sat, ...
                r_star, time_step)
            % Reduce data to valid occultation events (satellite and star
            % vectors must have angle > pi /2 between them
            dot_vec = dot(r_sat, r_star, 2);

            time = time(dot_vec < 0);

            if ~isempty(time)
                r_sat = r_sat(dot_vec < 0, :);
                r_star = r_star(dot_vec < 0, :);
    
                % Compute zenith angle (phi) and angle between satellite and
                % occultation point vectors (theta)
                cross_vec = cross(r_sat, r_star, 2);
                norm_cross_vec = vecnorm(cross_vec, 2, 2);
    
                phi = atan2(norm_cross_vec, dot_vec(dot_vec < 0));
                theta = phi - (pi / 2);
    
                % Calculate the occultation point altitude
                point_alt = vecnorm(r_sat, 1, 2) .* cos(theta);
    
                % Calculate the rotation axis
                r_hat = cross_vec ./ norm_cross_vec;
    
                % Normalize satellite vectors to be rotated
                v = r_sat ./ vecnorm(r_sat, 1, 2);
    
                % Calculate occultation point unit vector in ICRF frame
                % Determine rotation matrix for rotating through theta to
                % vector perpendicular to satellite-to-star vector. Adapted
                % from Analytical Mechanics of Space Systems Sec 3.4 by
                % Hanspeter Schaub and John L. Junkins
                p_hat = (((1 - cos(theta)) .* dot(v, r_hat, 2)) .* r_hat + ...
                    cos(theta) .* v + sin(theta) .* cross(r_hat, v, 2));
                point_vec = point_alt .* p_hat;
                
                % Calculate angle between occultation point and target star,
                % satellite and occultation point
                p_star_angle = acos(dot(p_hat, r_star ./ ...
                    vecnorm(r_star, 1, 2), 2));
                theta_act = atan2(vecnorm(cross(r_sat, point_vec, 2), 2, 2),...
                    dot(r_sat, point_vec, 2));
    
                % Find the break points in the time data to identify the
                % different occultation events (complete star-rises, star-sets)
                diff_check = logical(diff(time) > time_step);
                
                % Concatenating true to beginning of array so first access
                % counts as beginning of occultation event. Concatenating true
                % to end of array so an index is provided to use in calculating
                % the end of the final access. Essentially a false occultation
                % start for use in indexing.
                s_idx = find([true diff_check' true]');
                
                occ_events = Occ_Event.empty(length(s_idx) - 1, 0);
    
                for i = 1 : length(s_idx) - 1
                    idx = s_idx(i) : s_idx(i + 1) - 1;
                    
                    % Must break up LLA calculations to prevent memory
                    % shortage at desired step sizes (< 1 sec). Would likely
                    % break for long occultations
                    % Using eci2lla MATLAB function to convert ICRF coordinates
                    % calculated above to Latitude, Longitude, and Altitude. In
                    % MATLAB, ICRF is used as the basis for ECI.
                    % https://www.mathworks.com/help/aerotbx/ug/ecef2eci.html
                    % https://www.mathworks.com/help/aerotbx/ug/
                    % eci2lla.html#d123e48861
                    lla_point = eci2lla(point_vec(idx, :) .* (10^3), ...
                        obj.dt_helper.epsec_2_date_vec(time(idx)));
                    
                    % Convert altitude to km from m provided by eci2lla
                    lla_point(:, 3) = lla_point(:, 3) .* (10^-3);
    
                    % Create the occ_event object
                    occ_events(i) = Occ_Event(star_name, time(idx), ...
                        lla_point, point_vec(idx, :), p_star_angle(idx),...
                        theta(idx), theta_act(idx), r_sat(idx, :), ...
                        obj.scenario_start);
                end
            else
                occ_events = [];
            end
        end

        %% Function to truncate occultation data down to region of 
        % interest.
        % Inputs:
        % occultations  : Occ_Event vector, vector of Occ_Events with data
        %                 to truncate
        % lat_min       : Float vector, specification of the minimum
        %                 latitudes for the regions of interest, degrees
        % lat_max       : Float vector, specification of the maximum
        %                 latitudes for the regions of interest, degrees
        % long_min      : Float vector, specification of the minimum
        %                 longitudes for the regions of interest, degrees
        % long_max      : Float vector, specification of the maximum
        %                 longitudes for the regions of interest, degrees
        % alt_min       : Float vector, specification of the minimum
        %                 altitudes for the regions of interest, km
        % alt_max       : Float vector, specification of the maximum
        %                 altitudes for the regions of interest, km                         
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        % Outputs:
        % trunc_occ     : Occ_Event vector, Occ_Events with all
        %                 out-of-region information removed
        function trunc_occ = truncate(~, occultations, lat_min, lat_max, ...
                long_min, long_max, alt_min, alt_max, start_time, stop_time)
            for i = 1 : length(occultations)
                idx_save = zeros(occultations(i).num_data(), 1);

                for j = 1 : length(lat_min)
                    lat_lim = [lat_min(j) lat_max(j)];
                    long_lim = [long_min(j) long_max(j)];
                    alt_lim = [alt_min(j) alt_max(j)];
                    idx = occultations(i).search_region_time(lat_lim,...
                        long_lim, alt_lim, [start_time, stop_time]);
                    idx_save = idx_save | idx;
                end
                
                occultations(i).delete_data(~idx_save);
            end

            emp_occ_idx = 0 == occultations.num_data();
            occultations(emp_occ_idx) = [];
            trunc_occ = occultations;
        end
    end
end