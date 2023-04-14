%% Nicholas Jones - njones31@vt.edu
% Satellite_Collection
% The satellite collection class contains the list of satellite objects
% used in the analysis.

classdef Satellite_Collection < handle
    properties (SetAccess = private)
        % List of satellite objects
        satellite_list = Satellite.empty;
    end
    
    methods
        %% Constructor for Satellite_Collection object. All satellite
        % initial classical elements are referenced to STK's TrueOfDate
        % reference frame.
        % name_vector       : String array, contains satellite name(s)
        % data_matrix       : Float matrix, contains data for satellites.
        %                     Should be formatted by column as:
        %                       1. Initial semimajor axis - km
        %                       2. Initial eccentricity
        %                       3. Initial inclination - deg
        %                       4. Initial Longitude of Ascending Node -
        %                          deg
        %                       5. Initial argument of perigee - deg
        %                       6. Initial true anomaly - deg
        %                       7. Drag coefficient
        %                       8. Area to mass ratio - m^2 kg^-1
        % scenario_start    : String, formatted string representing the
        %                     scenario start time, follows format specified
        %                     in Datetime_Helper constructor.
        function obj = Satellite_Collection(name_vector, data_matrix, ...
                scenario_start)
            obj.add_sat(name_vector, data_matrix, scenario_start);
        end

        % Retrieves the satellite list from the collection
        % Outputs:
        % sat_list  : Satellite vector, vector of all satellites in the
        % collection
        function sat_list = get.satellite_list(obj)
            sat_list = obj.satellite_list;
        end

        %% Get the names of all satellites in the collection.
        % Outputs:
        % name_vector   : String array, contains all the satellite names
        function name_vector = get_name_vector(obj)
            name_vector = [obj.satellite_list.name];
        end

        %% Get the initial classical elements for all satellites in the
        % collection.
        % Outputs:
        % data_matrix   : Float vector, contains data for all satellites in
        %                 the satellite collection.
        function data_matrix = get_classical_elements(obj)
            data_matrix = [obj.satellite_list.a_init...
                obj.satellite_list.e_init obj.satellite_list.i_init ...
                obj.satellite_list.lan_init ...
                obj.satellite_list.omega_init obj.satellite_list.nu_init];
        end

        %% Get the aerodynamic properties of all satellites in the
        % collection.
        % Outputs:
        % data_matrix   : Float vector, contains drag coefficient and area
        %                 to mass ratio for all satellites in the
        %                 collection
        function data_matrix = get_aero_prop(obj)
            data_matrix = [obj.satellite_list.drag_coeff...
                obj.satellite_list.area_mass_ratio];
        end

        %% Add a satellite to the collection
        % Inputs:
        % name              : Strng array, name(s) of the satellite(s) to
        %                     add to the collection. If inputting a string
        %                     array, should be a column.
        % data              : Float matrix, satellite initial classical
        %                     elements formatted consistently with
        %                     constructor format
        % scenario_start    : String, formatted string representing the
        %                     scenario start time, follows format specified
        %                     in Datetime_Helper constructor.
        function add_sat(obj, name, data, scenario_start)
            if size(name, 1) ~= size(data, 1)
                error(['Satellite_Collection::add_sat: Length of names',...
                    ' does not match data']);
            elseif size(data, 2) ~= 8
                error(['Satellite_Collection::add_sat: Incorrect data'...
                    'formatting']);
            end

            for i = 1 : size(name, 1)
                if isempty(obj.satellite_list) ||...
                        ~any(matches([obj.satellite_list.name], name(i)))
                    obj.satellite_list(end + 1, 1) = ...
                        Satellite(scenario_start, name(i), data(i, 1), ...
                        data(i, 2), data(i, 3), data(i, 4), data(i, 5), ...
                        data(i, 6), data(i, 7), data(i, 8));
                else
                    error(['Satellite_Collection::add_sat: Duplicate', ...
                        ' name detected. Failed at: ' num2str(i)]);
                end
            end
        end

        %% Remove a satellite from the collection.
        % Inputs:
        % name  : String array, name(s) of the satellite(s) to remove from
        %         the collection. If inputting a string array, should be a
        %         column.
        function remove_sat(obj, name)
            remove_idx = matches([obj.satellite_list.name]', name);
            obj.satellite_list(remove_idx) = [];
        end

        %% Search the satellite collection by satellite names.
        % Inputs:
        % name      : String array, name(s) of the satellite(s) to search
        %             the collection for.
        % Outputs:
        % sat_list  : Satellite array, contains the found satellites.
        function [sat_list] = search_names(obj, name)
            search_idx = matches([obj.satellite_list.name], name);
            sat_list = obj.satellite_list(search_idx);
        end

        %% Search the satellite collection by initial semi-major axis in
        % km. Search limits are inclusive.
        % Inputs:
        % a_lim     : Float vector, specification of the min and max
        %             initial semimajor axis in km. [min_a max_a]
        % Outputs:
        % sat_list  : Satellite array, contains found satellites
        function sat_list = search_a(obj, a_lim)
            search_idx = [obj.satellite_list.a_init] >= a_lim(1) &...
                [obj.satellite_list.a_init] <= a_lim(2);
            sat_list = obj.satellite_list(search_idx);
        end

        %% Search the satellite collection by initial inclination in deg.
        % Search limits are inclusive.
        % Inputs:
        % i_lim     : Float vector, specification of the min and max
        %             initial inclination in deg. [min_i max_i]
        % Outputs:
        % sat_list  : Satellite array, contains found satellites
        function sat_list = search_i(obj, i_lim)
            search_idx = [obj.satellite_list.i_init] >= i_lim(1) &...
                [obj.satellite_list.i_init] <= i_lim(2);
            sat_list = obj.satellite_list(search_idx);
        end

        %% Search the satellite collection by initial Longitude of
        % Ascending Node in deg. Search limits are inclusive.
        % Inputs:
        % lan_lim   : Float vector, specification of the min and max
        %             initial Longitude of Ascending Node in deg.
        %             [min_lan max_alan]
        % Outputs:
        % sat_list  : Satellite array, contains found satellites
        function sat_list = search_lan(obj, lan_lim)
            search_idx = [obj.satellite_list.lan_init] >= lan_lim(1) &...
                [obj.satellite_list.lan_init] <= lan_lim(2);
            sat_list = obj.satellite_list(search_idx);
        end

        %% Function to get the number of satellites in the collection
        % Outputs:
        % num_sat   : Int, number of satellites in the collection
        function num_sat = count(obj)
            num_sat = length(obj.satellite_list);
        end

        %% Function to get a satellite at an index in the collection
        % Inputs:
        % idx   : Int, satellite index
        % Outputs:
        % sat   : Satellite, satellite at requested index
        function sat = get_sat(obj, idx)
            sat = obj.satellite_list(idx);
        end
    end
end