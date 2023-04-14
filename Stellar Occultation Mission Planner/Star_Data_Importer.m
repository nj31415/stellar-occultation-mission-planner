%% Nicholas Jones - njones31@vt.edu
% Star_Data_Importer
% Reads a csv file with star data and provides formatted matrices that can
% be supplied to Star_Collection to add Star objects.

classdef Star_Data_Importer < handle
    properties (SetAccess = public)
        data_file
    end
    
    methods
        %% Constructor for Star_Data_Importer class. Used to create a
        % Star_Data_Importer object to parse csv files.
        % Input:
        % data_file : String, file path to csv file containing the star
        %             data. File specified by data file should have the
        %             following format by column:
        %               1. Star name
        %               2. HD catalog number
        %               3. Right ascension - degrees
        %               4. Declination - degrees
        %               5. Right ascension proper motion - arcsec yr^-1
        %               6. Declination proper motion - arcsec yr^-1
        %               7. Parallax - arcsec
        function obj = Star_Data_Importer(data_file)
            obj.data_file = data_file;
        end
        
        %% Get the data file path being used by Star_Data_Importer
        % Outputs:
        % data_file : String, file path to csv file containing the star
        %             data.
        function data_file = get.data_file(obj)
            data_file = obj.data_file;
        end

        %% Set the data file path being used by Star_Data_Importer.
        % Inputs:
        % data_file : String, file path to csv file containing the star
        %             data.
        function set.data_file(obj, data_file)
            if isfile(data_file)
                obj.data_file = data_file;
            else
                error(['Star_Data_Importer::set.data_file: Could not ', ...
                    'find the specified file.']);
            end
        end

        %% Parse the data file and extract star data.
        % Outputs:
        % name_vector   : String array, containing star name(s).
        % data_matrix   : Float matrix, contains data for stars.
        function [name_vector, data_matrix] = parse_file(obj)
            data = readtable(obj.data_file, 'VariableNamingRule', 'preserve');
            name_vector = table2array(data(:, 1));
            
            % Replace spaces with underscore in identifiers
            name_vector = strrep(name_vector, ' ', '_');
            data_matrix = table2array(data(:, 2 : end));
        end
    end
end