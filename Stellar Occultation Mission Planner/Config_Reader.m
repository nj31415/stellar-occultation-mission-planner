%% Nicholas Jones - njones31@vt.edu
% The Config_Reader class is used to read in configuration options from a
% config file.

classdef Config_Reader < handle
    properties (SetAccess = immutable)
        config_file
    end
    
    methods
        %% Constructor for a Config_Reader object. Specify the config file
        % to extract options from.
        % Inputs:
        % config_file   : String, filepath to config text file.
        function obj = Config_Reader(config_file)
            if isfile(config_file)
                obj.config_file = config_file;
            else
                error(['Config_Reader::Constructor: Could not find ',...
                    'indicated file']);
            end
        end
        
        %% Get the file path of the config file being used for the
        % Config_Reader.
        % Outputs:
        % config_file   : String, config file filepath
        function config_file = get.config_file(obj)
            config_file = obj.config_file;
        end

        %% Parse config file
        % Outputs:
        % save_folder           : String, filepath to folder to save data
        %                         to.
        % run_options_folder    : String, filepath to folder containing run
        %                         options for the simulation.
        % star_data_folder      : String, filepath to folder containing
        %                         star data.
        % sat_data_folder       : String, filepath to folder containing
        %                         satellite data
        % stk_version           : int, STK version number
        % stk_scenario_folder   : String, filepath to folder containing the
        %                         STK scenario files.
        % stk_visible           : boolean, controls visibility of STK
        % run_select            : String, specifies run_options file to use
        %                         for the simulation.
        function [save_folder, run_options_folder, star_data_folder, ...
                sat_data_folder, stk_version, stk_scenario_folder, ...
                stk_visible, run_select] = parse_file(obj)
            lines = readlines(obj.config_file);

            % Remove comment lines (indicated by leading '%' character).
            idx = [];
            for i = 1 : length(lines)
                if isempty(lines{i}) || lines{i}(1) == '%'
                    idx(end + 1) = i;
                end
            end
            lines(idx) = [];

            save_folder = strsplit(lines{1}, ' = ');
            save_folder = save_folder{2}(1 : end - 1);

            run_options_folder = strsplit(lines{2}, ' = ');
            run_options_folder = run_options_folder{2}(1 : end - 1);

            star_data_folder = strsplit(lines{3}, ' = ');
            star_data_folder = star_data_folder{2}(1 : end - 1);

            sat_data_folder = strsplit(lines{4}, ' = ');
            sat_data_folder = sat_data_folder{2}(1 : end - 1);

            stk_version = strsplit(lines{5}, ' = ');
            stk_version = str2double(stk_version{2}(1 : end - 1));

            stk_scenario_folder = strsplit(lines{6}, ' = ');
            stk_scenario_folder = stk_scenario_folder{2}(1 : end - 1);

            stk_visible = strsplit(lines{7}, ' = ');
            stk_visible = str2num(stk_visible{2}(1 : end - 1));

            run_select = strsplit(lines{8}, ' = ');
            run_select = run_select{2}(1 : end - 1);
        end
    end
end