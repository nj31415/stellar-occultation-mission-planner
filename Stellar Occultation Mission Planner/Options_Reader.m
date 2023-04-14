%% Nicholas Jones - njones31@vt.edu
% The Options_Reader class is used to read in simulation options for a
% specific run.

classdef Options_Reader < handle
    properties (SetAccess = immutable)
        options_file
    end
    
    methods
        %% Constructor for an Options_Reader object. Specify the options
        % file to extract options from.
        % Inputs:
        % run_options   : String, filepath to options text file.
        function obj = Options_Reader(run_options)
            if isfile(run_options)
                obj.options_file = run_options;
            else
                error(['Options_Reader::Constructor: Could not find ',...
                    'indicated file']);
            end
        end
        
        %% Get the file path of the options file being used for the
        % Options_Reader.
        % Outputs:
        % options_file   : String, config file filepath
        function options_file = get.options_file(obj)
            options_file = obj.options_file;
        end

        %% Parse the options file
        % Outputs:
        % run_name          : String, name of the run. Will be used in save
        %                     file generation if enabled.
        % load_data         : Boolean, controls whether to load data or
        %                     simulate using STK
        % save_data         : Boolean, controls whether to save data
        % import_stars      : Boolean, controls whether to import stars
        % import_sats       : Boolean, controls whether to import
        %                     satellites
        % truncate_data     : Boolean, controls whether to truncate data to
        %                     region of interest while saving.
        % plot_data         : Boolean, controls whether to plot data at the
        %                     end of the run.
        % stk_scenario_file : String, STK scenario file to use.
        % star_data_file    : String, star data file to use.
        % sat_data_file     : String, satellite data file to use.
        % load_data_file    : String, save data file to use if not
        %                     simulating in STK.
        % ephem_step_size   : Float, ephemeris extraction step size in sec
        % occ_step_size     : Float, access calculation step size in sec
        % start_time        : String, formatted date time to set the
        %                     scenario start time. Follows format:
        %                     yyyy:mm:dd:HH:MM:SS.FFF
        % stop_time         : String, formatted date time to set the
        %                     scenario end time. Follows same format as
        %                     start_time.
        % oa_start_time     : String, formatted date time to set the
        %                     analysis start time for computing
        %                     occultations. Follows same format as
        %                     start_time, may be empty.
        % oa_stop_time      : String, formatted date time to set the
        %                     analysis stop time for computing
        %                     occultations. Follows same format as
        %                     start_time, may be empty. 
        % propagator        : String, propagator to use in STK.
        % sat_names         : String array, name of satellites to perform
        %                     the analysis for. Should match entries in
        %                     satellite data file. If empty, use all
        %                     satellites in the satellite data file.
        % star_names        : String array, name of the stars to calculate
        %                     occultations for. Should match entries in
        %                     star data file. If empty, use all stars in
        %                     the star data file.
        % lat_min           : Float vector, minimum latitudes of region of
        %                     interest. Multiple entries indicate multiple
        %                     regions of interest, aligned by index.
        % lat_max           : Float vector, maximum latitudes of region of
        %                     interest. Multiple entries indicate multiple
        %                     regions of interest, aligned by index.
        % long_min          : Float vector, minimum longitudes of region of
        %                     interest. Multiple entries indicate multiple
        %                     regions of interest, aligned by index.
        % long_max          : Float vector, maximum longitudes of region of
        %                     interest. Multiple entries indicate multiple
        %                     regions of interest, aligned by index.
        % alt_min           : Float vector, minimum altitudes of region of
        %                     interest. Multiple entries indicate multiple
        %                     regions of interest, aligned by index.
        % alt_max           : Float vector, maximum altitudes of region of
        %                     interest. Multiple entries indicate multiple
        %                     regions of interest, aligned by index.
        % max_star_angle    : Float, maximum angle between the star and
        %                     orbit plane to be considered a valid
        %                     occultation.
        % min_sun_angle     : Float, minimum angle between the instrument
        %                     (vector from spacecraft to star) and the Sun
        %                     to be considered a valid occultation.
        % min_moon_angle    : Float, minimum angle between the instrument
        %                     (vector from spacecraft to star) and the Moon
        %                     to be considered a valid occultation.
        % atmos_model       : Formatted string for setting the HPOP
        %                     atmospheric model
        function [run_name, load_data, save_data, import_stars, ...
                import_sats, truncate_data, plot_data, ...
                stk_scenario_file, star_data_file, sat_data_file, ...
                load_data_file, ephem_step_size, occ_step_size, ...
                start_time, stop_time, oa_start_time, oa_stop_time, ...
                propagator, sat_names, star_names, lat_min, lat_max,...
                long_min, long_max, alt_min, alt_max, max_star_angle, ...
                min_sun_angle, min_moon_angle, atmos_model] = ...
                parse_file(obj)
            lines = readlines(obj.options_file);

            % Remove comment lines (indicated by leading '%' character).
            idx = [];
            for i = 1 : length(lines)
                if isempty(lines{i}) || lines{i}(1) == '%'
                    idx(end + 1) = i;
                end
            end
            lines(idx) = [];

            run_name = strsplit(lines{1}, ' = ');
            run_name = run_name{2}(1 : end - 1);

            load_data = strsplit(lines{2}, ' = ');
            load_data = str2num(load_data{2}(1 : end - 1));

            save_data = strsplit(lines{3}, ' = ');
            save_data = str2num(save_data{2}(1 : end - 1));

            import_stars = strsplit(lines{4}, ' = ');
            import_stars = str2num(import_stars{2}(1 : end - 1));

            import_sats = strsplit(lines{5}, ' = ');
            import_sats = str2num(import_sats{2}(1 : end - 1));

            truncate_data = strsplit(lines{6}, ' = ');
            truncate_data = str2num(truncate_data{2}(1 : end - 1));

            plot_data = strsplit(lines{7}, ' = ');
            plot_data = str2num(plot_data{2}(1 : end - 1));

            stk_scenario_file = strsplit(lines{8}, ' = ');
            stk_scenario_file = stk_scenario_file{2}(1 : end - 1);

            star_data_file = strsplit(lines{9}, ' = ');
            star_data_file = star_data_file{2}(1 : end - 1);

            sat_data_file = strsplit(lines{10}, ' = ');
            sat_data_file = sat_data_file{2}(1 : end - 1);

            load_data_file = strsplit(lines{11}, ' = ');
            load_data_file = load_data_file{2}(1 : end - 1);

            ephem_step_size = strsplit(lines{12}, ' = ');
            ephem_step_size = str2double(ephem_step_size{2}(1 : end - 1));

            occ_step_size = strsplit(lines{13}, ' = ');
            occ_step_size = str2double(occ_step_size{2}(1 : end - 1));

            start_time = strsplit(lines{14}, ' = ');
            start_time = start_time{2}(1 : end - 1);

            stop_time = strsplit(lines{15}, ' = ');
            stop_time = stop_time{2}(1 : end - 1);

            oa_start_time = strsplit(lines{16}, ' = ');
            oa_start_time = oa_start_time{2}(1 : end - 1);

            oa_stop_time = strsplit(lines{17}, ' = ');
            oa_stop_time = oa_stop_time{2}(1 : end - 1);

            propagator = strsplit(lines{18}, ' = ');
            propagator = propagator{2}(1 : end - 1);

            sat_names = strsplit(lines{19}, ' = ');
            sat_names = strsplit(sat_names{2}(1 : end - 1), ', ');

            if isempty(sat_names{1})
                sat_names = [];
            end

            star_names = strsplit(lines{20}, ' = ');
            star_names = strsplit(star_names{2}(1 : end - 1), ', ');

            if isempty(star_names{1})
                star_names = [];
            end

            lat_min = strsplit(lines{21}, ' = ');
            lat_min = str2double(strsplit(lat_min{2}(1 : end - 1), ', '));

            lat_max = strsplit(lines{22}, ' = ');
            lat_max = str2double(strsplit(lat_max{2}(1 : end - 1), ', '));

            long_min = strsplit(lines{23}, ' = ');
            long_min = str2double(strsplit(long_min{2}(1 : end - 1), ...
                ', '));

            long_max = strsplit(lines{24}, ' = ');
            long_max = str2double(strsplit(long_max{2}(1 : end - 1), ...
                ', '));

            alt_min = strsplit(lines{25}, ' = ');
            alt_min = str2double(strsplit(alt_min{2}(1 : end - 1), ', '));

            alt_max = strsplit(lines{26}, ' = ');
            alt_max = str2double(strsplit(alt_max{2}(1 : end - 1), ', '));

            max_star_angle = strsplit(lines{27}, ' = ');
            max_star_angle = str2double(max_star_angle{2}(1 : end - 1));

            min_sun_angle = strsplit(lines{28}, ' = ');
            min_sun_angle = str2double(min_sun_angle{2}(1 : end - 1));

            min_moon_angle = strsplit(lines{29}, ' = ');
            min_moon_angle = str2double(min_moon_angle{2}(1 : end - 1));
            
            atmos_model = strsplit(lines{30}, ' = ');
            atmos_model = atmos_model{2}(1 : end - 1);
        end
    end
end