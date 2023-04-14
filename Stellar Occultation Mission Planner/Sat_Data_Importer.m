%% Nicholas Jones - njones31@vt.edu
% Sat_Data_Importer
% Reads a csv file with satellite data and provides formatted matrices that
% can be supplied to Satellite_Collection to add Satellite objects.

classdef Sat_Data_Importer < handle
    properties (SetAccess = public)
        data_file
    end
    
    methods
        %% Constructor for Sat_Data_Importer class. Used to create a
        % Sat_Data_Importer object to parse csv files.
        % Input:
        % data_file : String, file path to csv file containing the
        %             satellite data. File specified by data file should
        %             have the following format by column:
        %               1. Satellite name
        %               2. Initial semi-major axis - km
        %               3. Initial eccentricity
        %               4. Initial inclination - deg
        %               5. Initial Longitude of Ascending Node - deg
        %               6. Initial argument of perigee - deg
        %               7. Initial true anomaly - deg
        %               8. Drag coefficient
        %               9. Area to mass ratio - m^2 kg^-1
        function obj = Sat_Data_Importer(data_file)
            obj.data_file = data_file;
        end
        
        %% Get the data file path being used by Sat_Data_Importer
        % Outputs:
        % data_file : String, file path to csv file containing the
        %             satellite data.
        function data_file = get.data_file(obj)
            data_file = obj.data_file;
        end

        %% Set the data file path being used by Star_Data_Importer.
        % Inputs:
        % data_file : String, file path to csv file containing the
        %             satellite data.
        function set.data_file(obj, data_file)
            if isfile(data_file)
                obj.data_file = data_file;
            else
                error(['Sat_Data_Importer::set.data_file: Could not ', ...
                    'find the specified file.']);
            end
        end

        %% Parse the data file and extract satellite data.
        % Outputs:
        % name_vector   : String array, containing satellite name(s).
        % data_matrix   : Float matrix, contains data for satellites.
        function [name_vector, data_matrix] = parse_file(obj)
            data = readtable(obj.data_file, 'VariableNamingRule', 'preserve');
            name_vector = table2array(data(:, 1));

            % Replace spaces with underscore in identifiers
            name_vector = strrep(name_vector, ' ', '_');
            data_matrix = table2array(data(:, 2 : end));
        end
    end
end