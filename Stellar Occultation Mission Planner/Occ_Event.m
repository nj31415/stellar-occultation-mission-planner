%% Nicholas Jones - njones31@vt.edu
% Occ_Event
% The Occ_Event class contains the occultation data between a star and a
% satellite. Each Occ_Event is between a specific star and satellite for a
% specific occultation event (complete star-rise or star-set). The
% Occ_Event class supports data retrieval based on the type of data
% requested, by time period, and by region of interest. Latitude, 
% longitude, and altitude of Occ_Events are referenced in the WGS-84
% system.

classdef Occ_Event < handle
    properties (SetAccess = public)
        star_name
    end

    properties (SetAccess = private)
        % data_matrix structure
        % 1.    Time Stamp - EpSec (seconds since scenario start time)
        % 2.    Latitude - deg
        % 3.    Longitude - deg
        % 4.    Altitude - km (from Earth's surface)
        % 5.    ICRF X - km
        % 6.    ICRF Y - km
        % 7.    ICRF Z - km
        % 8.    Angle between tangent point position and star vector - rad
        %       (should be pi / 2 rad)
        % 9.    Required angle between tangent point position and satellite
        %       position vector (theta) - rad
        % 10.   Actual angle between point position and satellite position
        %       vector (theta) - rad
        % 11.   Satellite ICRF X - km
        % 12.   Satellite ICRF Y - km
        % 13.   Satellite ICRF Z - km
        data_matrix = zeros(1, 13);

        dt_helper   % Datetime_Helper object
    end
    
    methods
        %% Constructor for an Occ_Event object. All geographic data is
        % referenced to WGS-84. All cartesian data is referenced to ICRF.
        % The properties of an Occ_Event can only be set through the
        % constructor, otherwise values are read-only.
        % Inputs:
        % star_name                 : String, the name of the star the
        %                             occultation data corresponds to.
        % time_data                 : Float vector, time stamps for the
        %                             occultation data in EpSec, referenced
        %                             to the beginning of the STK scenario
        %                             period.
        % lla_data                  : Float matrix, latitude, longitude,
        %                             and altitude in deg and km
        %                             respectively.
        % icrf_data                 : Float matrix, ICRF X, Y, and Z
        %                             position data in km.
        % angle_pos_star_data       : Float vector, angle between the
        %                             tangent point position and star
        %                             vector in rad. Should be pi / 2 rad.
        % angle_pos_sat_data_req    : Float vector, angle between the
        %                             tangent point position and satellite 
        %                             position from theory in rad.
        % angle_pos_sat_data_act    : Float vector, angle between the
        %                             tangent point position and satellite 
        %                             position as implemented in the model
        %                             in rad.
        % sat_icrf_data             : Float matrix, ICRF X, Y, and Z
        %                             position of the satellite in km
        % scenario_start            : String, formatted string representing
        %                             the scenario start time, follows
        %                             format specified in Datetime_Helper
        %                             constructor.
        function obj = Occ_Event(star_name, time_data, lla_data, ...
                icrf_data, angle_pos_star_data, angle_pos_sat_data_req, ...
                angle_pos_sat_data_act, sat_icrf_data, scenario_start)
            obj.star_name = star_name;
            
            n = length(time_data);

            obj.data_matrix(1 : n, 1) = time_data;
            obj.data_matrix(1 : n, 2 : 4) = lla_data;
            obj.data_matrix(1 : n, 5 : 7) = icrf_data;
            obj.data_matrix(1 : n, 8) = angle_pos_star_data;
            obj.data_matrix(1 : n, 9) = angle_pos_sat_data_req;
            obj.data_matrix(1 : n, 10) = angle_pos_sat_data_act;
            obj.data_matrix(1 : n, 11 : 13) = sat_icrf_data;

            obj.dt_helper = Datetime_Helper(scenario_start);
        end
        
        %% Get the name of the star corresponding to this Occ_Event.
        % Outputs:
        % star_name : String, the name of the star corresponding to the
        %             Occ_event
        function star_name = get.star_name(obj)
            star_name = obj.star_name;
        end
        
        %% Set the name of the star corresponding to this Occ_Event. This
        % must be the name of a star loaded into the Star_Collection for
        % the simulation.
        % Inputs:
        % star_name : String, the star name
        function set.star_name(obj, star_name)
            obj.star_name = star_name;
        end

        %% Get the time stamp data for the Occ_Event.
        % Inputs:
        % ind       : Index vector, the desired indices to retrieve. If
        %             unspecified, retrieves all data.
        % Outputs:
        % time_data : Float vector, time stamp data for the occultation
        %             between the satellite and star in EpSec.
        function time_data = get_time_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                time_data = obj.data_matrix(ind, 1);
            else
                time_data = obj.data_matrix(:, 1);
            end
        end

        %% Get the latitude data for the Occ_Event.
        % Inputs:
        % ind       : Index vector, the desired indices to retrieve. If
        %             unspecified, retrieves all data.
        % Outputs:
        % lat_data  : Float vector, latitude data for the occultation
        %             between the satellite and star in degrees.
        function lat_data = get_lat_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                lat_data = obj.data_matrix(ind, 2);
            else
                lat_data = obj.data_matrix(:, 2);
            end
        end

        %% Get the longitude data for the Occ_Event.
        % Inputs:
        % ind       : Index vector, the desired indices to retrieve. If
        %             unspecified, retrieves all data.
        % Outputs:
        % long_data : Float vector, longitude data for the occultation
        %             between the satellite and star in degrees.
        function long_data = get_long_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                long_data = obj.data_matrix(ind, 3);
            else
                long_data = obj.data_matrix(:, 3);
            end
        end

        %% Get the altitude data for the Occ_Event.
        % Inputs:
        % ind       : Index vector, the desired indices to retrieve. If
        %             unspecified, retrieves all data.
        % Ouputs:
        % alt_data  : Float vector, altitude data for the occultation
        %             between the satellite and star in km above the
        %             Earth's surface.
        function alt_data = get_alt_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                alt_data = obj.data_matrix(ind, 4);
            else
                alt_data = obj.data_matrix(:, 4);
            end
        end

        %% Get the LLA and time stamp data for the Occ_Event.
        % Inputs:
        % ind       : Index vector, the desired indices to retrieve. If
        %             unspecified, retrieves all data.
        % Ouputs:
        % lla_data  : Float vector, the time stamp, latitude, longitude,
        %             and altitude data for the occultation between the
        %             star in a single matrix.
        function lla_data = get_lla_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                lla_data = obj.data_matrix(ind, 1 : 4);
            else
                lla_data = obj.data_matrix(:, 1 : 4);
            end
        end

        %% Get the ICRF X position for the Occ_Event.
        % Inputs:
        % ind           : Index vector, the desired indices to retrieve. If
        %                 unspecified, retrieves all data.
        % Ouputs:
        % icrf_x_data   : Float vector, the ICRF X position in km for the
        %                 occultation between the satellite and the star
        function icrf_x_data = get_icrf_x_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                icrf_x_data = obj.data_matrix(ind, 5);
            else
                icrf_x_data = obj.data_matrix(:, 5);
            end
        end

        %% Get the ICRF Y position for the Occ_Event.
        % Inputs:
        % ind           : Index vector, the desired indices to retrieve. If
        %                 unspecified, retrieves all data.
        % Ouputs:
        % icrf_y_data   : Float vector, the ICRF Y position in km for the
        %                 occultation between the satellite and the star
        function icrf_y_data = get_icrf_y_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                icrf_y_data = obj.data_matrix(ind, 6);
            else
                icrf_y_data = obj.data_matrix(:, 6);
            end
        end

        %% Get the ICRF Z position for the Occ_Event.
        % Inputs:
        % ind           : Index vector, the desired indices to retrieve. If
        %                 unspecified, retrieves all data.
        % Ouputs:
        % icrf_z_data   : Float vector, the ICRF Z position in km for the
        %                 occultation between the satellite and the star
        function icrf_z_data = get_icrf_z_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                icrf_z_data = obj.data_matrix(ind, 7);
            else
                icrf_z_data = obj.data_matrix(:, 7);
            end
        end

        %% Get the ICRF Position and time stamp data for the Occ_Event.
        % Inputs:
        % ind       : Index vector, the desired indices to retrieve. If
        %             unspecified, retrieves all data.
        % Ouputs:
        % icrf_data : Float vector, the time stamp and ICRF position data
        %             for the occultation between the star in a single
        %             matrix.
        function icrf_data = get_icrf_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                icrf_data = obj.data_matrix(ind, [1, 5 : 7]);
            else
                icrf_data = obj.data_matrix(:, [1, 5 : 7]);
            end
        end

        %% Get the satellite ICRF X position for the Occ_Event.
        % Inputs:
        % ind           : Index vector, the desired indices to retrieve. If
        %                 unspecified, retrieves all data.
        % Ouputs:
        % icrf_x_data   : Float vector, the ICRF X position in km for the
        %                 occultation between the satellite and the star
        function sat_icrf_x_data = get_sat_icrf_x_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                sat_icrf_x_data = obj.data_matrix(ind, 11);
            else
                sat_icrf_x_data = obj.data_matrix(:, 11);
            end
        end

        %% Get the satellite ICRF Y position for the Occ_Event.
        % Inputs:
        % ind           : Index vector, the desired indices to retrieve. If
        %                 unspecified, retrieves all data.
        % Ouputs:
        % icrf_y_data   : Float vector, the ICRF Y position in km for the
        %                 occultation between the satellite and the star
        function sat_icrf_y_data = get_sat_icrf_y_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                sat_icrf_y_data = obj.data_matrix(ind, 12);
            else
                sat_icrf_y_data = obj.data_matrix(:, 12);
            end
        end

        %% Get the satellite ICRF Z position for the Occ_Event.
        % Inputs:
        % ind           : Index vector, the desired indices to retrieve. If
        %                 unspecified, retrieves all data.
        % Ouputs:
        % icrf_z_data   : Float vector, the ICRF Z position in km for the
        %                 occultation between the satellite and the star
        function sat_icrf_z_data = get_sat_icrf_z_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                sat_icrf_z_data = obj.data_matrix(ind, 13);
            else
                sat_icrf_z_data = obj.data_matrix(:, 13);
            end
        end

        %% Get the satellite ICRF Position and time stamp data for the
        % Occ_Event.
        % Inputs:
        % ind       : Index vector, the desired indices to retrieve. If
        %             unspecified, retrieves all data.
        % Ouputs:
        % icrf_data : Float vector, the time stamp and ICRF position data
        %             for the occultation between the star in a single
        %             matrix.
        function sat_icrf_data = get_sat_icrf_data(obj, ind)
            if exist('ind', 'var') && ~isempty(ind)
                sat_icrf_data = obj.data_matrix(ind, [1, 11 : 13]);
            else
                sat_icrf_data = obj.data_matrix(:, [1, 11 : 13]);
            end
        end

        %% Function to delete data at provided indices. Should check if all
        % data from occultation is deleted, and if so, delete the Occ_Event
        % Inputs:
        % idx   : index vector, deletes all data at specified indices.
        function delete_data(obj, idx)
            obj.data_matrix(idx, :) = [];
        end

        %% Function to retrieve length of data matrix for Occ_Event object
        % Outputs:
        % length    : Int, number of entries in data matrix
        function num = num_data(obj)
            num = zeros(size(obj, 2), 1);
            for i = 1 : length(obj)
                num(i) = size(obj(i).data_matrix, 1);
            end
        end

        %% Retrieve logical indices for tangent points within a given
        % geographic region. Searches are inclusive of the limits. 
        % Can be used as an input for the getter functions within this
        % class.
        % Inputs:
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % Outputs:
        % ind_region    : Logical index array, true where the tangent point
        %                 meets the region of interest conditions.
        function ind_region = ...
                search_region(obj, lat_lim, long_lim, alt_lim)
            ind_region = obj.data_matrix(:, 2) >= lat_lim(1) & ...
                obj.data_matrix(:, 2) <= lat_lim(2) & ...
                obj.data_matrix(:, 3) >= long_lim(1) & ...
                obj.data_matrix(:, 3) <= long_lim(2) & ...
                obj.data_matrix(:, 4) >= alt_lim(1) & ...
                obj.data_matrix(:, 4) <= alt_lim(2);
        end

        %% Retrieve logical indices for tangent points within a given time
        % interval. Searches are inclusive of the limits. Can be used as an
        % input for the getter functions within this class.
        % Inputs:
        % time_lim  : Float vector, specification of the minimum and
        %             maximum times in EpSec. [min_time max_time]
        % Outputs:
        % ind_time  : Logical index array, true when the tangent point
        %             occurs within the time period of interest.
        function ind_time = search_time(obj, time_lim)
            ind_time = obj.data_matrix(:, 1) >= time_lim(1) & ...
                obj.data_matrix(:, 1) <= time_lim(2);
        end

        %% Retrieve logical indices for tangent points within a given time 
        % interval and within a geographic region of interest. Searches are
        % inclusive of the limits. Can be used as an input for the getter
        % functions within this class.
        % Inputs:
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % time_lim      : Float vector, specification of the minimum and
        %                 maximum times in EpSec. [min_time max_time]
        % Outputs:
        % ind           : Logical index array, true where the tangent point
        %                 meets the region of interest and time conditions.
        function ind = search_region_time(obj, ...
                lat_lim, long_lim, alt_lim, time_lim)
            ind = obj.data_matrix(:, 2) >= lat_lim(1) & ...
                obj.data_matrix(:, 2) <= lat_lim(2) & ...
                obj.data_matrix(:, 3) >= long_lim(1) & ...
                obj.data_matrix(:, 3) <= long_lim(2) & ...
                obj.data_matrix(:, 4) >= alt_lim(1) & ...
                obj.data_matrix(:, 4) <= alt_lim(2) & ...
                obj.data_matrix(:, 1) >= time_lim(1) & ...
                obj.data_matrix(:, 1) <= time_lim(2);
        end
        
        %% Retrieve the start time of the occultation.
        % Outputs:
        % time_start    : Float, start time of the occultation in EpSec
        function time_start = get_start_time(obj)
            time_start = obj.data_matrix(1, 1);
        end

        %% Retrieve the end time of the occultation.
        % Outputs:
        % time_end      : Float, end time of the occultation in EpSec
        function time_end = get_time_end(obj)
            time_end = obj.data_matrix(end, 1);
        end

        %% Retrieve the ICRF velocity vector for the occultation event
        % tangent points in km / s. Calculated by forward differencing of
        % the ICRF position data. If the last point is the final point in
        % the occultation event, velocity will be calculated using backward
        % differencing.
        % Inputs:
        % ind           : Index vector, the desired indices to retrieve. If
        %                 unspecified, retrieves all data.
        % Ouputs:
        % vel_icrf      : Float matrix, ICRF velocites in km / s. First
        %                 column is x, then y, then z axes.
        function vel_icrf = get_icfr_vel(obj, ind)
            if any(ind)
                icrf_data = obj.get_icrf_data(ind);
                icrf_dif = diff(icrf_data, 1, 1);
                vel_icrf = icrf_dif(:, 2 : 4) ./ icrf_dif(:, 1);
                
                if ~isempty(vel_icrf)
                    % Calculate the velocity of the last data point
                    if islogical(ind)
                        ind_l_p1 = find(ind, 1, 'last') + 1;
                    else
                        ind_l_p1 = ind(end) + 1;
                    end
        
                    if ind_l_p1 <= size(obj.data_matrix, 1)
                        % Forward differencing calculation for the last
                        % desired index, using data from the next index if 
                        % it exists
                        icrf_dif_l = (obj.get_icrf_data(ind_l_p1) - ...
                            icrf_data(end, :));
                        vel_icrf(:, end + 1) = icrf_dif_l(2 : 4) ./ ...
                            icrf_dif_l(1);
                    else
                        % Backward differencing is the same as forward 
                        % differencing from the previous index
                        vel_icrf(:, end + 1) = vel_icrf(:, end);
                    end
                else
                    vel_icrf = zeros(3, 0);
                end
            else
                vel_icrf = double.empty(3, 0);
            end
        end

        %% Retrieve the altitude rate of change for the occultation event
        % tangent points in km / s. Calculated by forward differencing of
        % the altitude LLA data. Returns empty if idx is emptys
        % Inputs:
        % ind       : Index vector, the desired indices to retrieve. If
        %             unspecified, retrieves all data.
        % Outputs:
        % alt_roc   : Float vector, rate of change of the altitude of the
        %             tangent point
        function alt_roc = get_alt_roc(obj, ind)
            if any(ind)
                alt_data = obj.get_alt_data(ind);
                time_data = obj.get_time_data(ind);
                alt_dif = diff([time_data alt_data], 1, 1);
                alt_roc = alt_dif(:, 2) ./ alt_dif(:, 1);
                if ~isempty(alt_roc)
                    % Calculate the velocity of the last data point
                    if islogical(ind)
                        ind_l_p1 = find(ind, 1, 'last') + 1;
                    else
                        ind_l_p1 = ind(end) + 1;
                    end
        
                    if ind_l_p1 <= size(obj.data_matrix, 1)
                        % Forward differencing calculation for the last
                        % desired index, using data from the next index if 
                        % it exists
                        alt_dif_l = (obj.get_alt_data(ind_l_p1) - ...
                            alt_data(end));
                        time_dif_l = (obj.get_time_data(ind_l_p1) - ...
                            time_data(end));
                        alt_roc(end + 1) = alt_dif_l / time_dif_l;
                    else
                        % Backward differencing is the same as forward
                        % differencing from the previous index
                        alt_roc(end + 1) = alt_roc(end);
                    end
                else
                    % Handle case where only 1 index is available and
                    % alt_roc cannot be calculated
                    alt_roc = 0;
                end

                alt_roc = reshape(alt_roc, [], 1);
            else
                alt_roc = [];
            end
        end

        %% Function to calculate the range over the duration of the
        % occultation between the satellite and occultation tangent point
        % Input:
        % ind       : Index vector, the desired indices to retrieve. If
        %             unspecified, retrieves all data.
        % Output:
        % range     : Float vector, range in km between the satellite and
        %             occultation tangent point
        function range = range_sat_tan(obj, ind)
            range = sqrt(...
                (obj.get_icrf_x_data(ind) - ...
                obj.get_sat_icrf_x_data(ind)).^2 + ...
                (obj.get_icrf_z_data(ind) - ...
                obj.get_sat_icrf_z_data(ind)).^2 + ...
                (obj.get_icrf_z_data(ind) - ...
                obj.get_sat_icrf_z_data(ind)).^2);
        end

        %% Function to calculate the mean local solar time of the
        % occultation tangent points
        function lst = local_solar_time(obj, ind)
            dt_data = ...
                obj.dt_helper.epsec_2_time(obj.get_time_data(ind));
            ut1_data = ...
                dt_data + seconds(deltaUT1(mjuliandate(dt_data)));

            % Adjust longitude data to between 0 and 360 degrees instead of
            % -180 to 180 degrees.
            long_data_2 = obj.get_long_data(ind);
            long_data_2(long_data_2 < 0) = long_data_2(long_data_2 < 0) ...
                + 360;

            lst = mod(hour(ut1_data) + (minute(ut1_data) ./ 60) + ...
                (second(ut1_data) ./ 3600) + (long_data_2 / 15), 24);
        end
    end
end