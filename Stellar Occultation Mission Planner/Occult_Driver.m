%% Nicholas Jones - njones31@vt.edu
% Driver script for a stellar occultation study.
close all;
clear;
clc;

%% Specify config file to use.
config = 'config.txt';

% Extract the config file options
cfg_reader = Config_Reader(config);
[save_folder, run_options_folder, star_data_folder, sat_data_folder,...
    stk_version, stk_scenario_folder, stk_visible, run_select] = ...
    cfg_reader.parse_file();

%% Extract the run options
opt_reader = Options_Reader([run_options_folder '/' run_select]);
[run_name, load_data, save_data, import_stars, import_sats, ...
    truncate_data, plot_data, stk_scenario_file, star_data_file, ...
    sat_data_file, load_data_file, ephem_step_size, occ_step_size, ...
    start_time, stop_time, oa_start_time, oa_stop_time, propagator, ...
    sat_names, star_names, lat_min, lat_max, long_min, long_max, ...
    alt_min, alt_max, max_star_angle, min_sun_angle, min_moon_angle, ...
    atmos_model] = opt_reader.parse_file();

% If the occultation analysis start and stop times are empty, sets them
% equal to the propagation start and stop times
if isempty(oa_start_time)
    oa_start_time = start_time;
end

if isempty(oa_stop_time)
    oa_stop_time = stop_time;
end

%% if load_data is true, retrieve simulation data from previous save file
% indicated in run options. Else, run simulation using STK
if load_data
    load([save_folder '/' load_data_file]);
else
    dt_helper = Datetime_Helper(start_time);

    tic;

    % Import stars into simulation
    star_importer = Star_Data_Importer(...
        [star_data_folder '/' star_data_file]);
    [star_name_vector, star_data_matrix] = star_importer.parse_file();
    star_collect = Star_Collection(star_name_vector, star_data_matrix);

    % If a list of stars had been specified in the run options, remove
    % those stars not on the list.
    if ~isempty(star_names)
        star_collect.remove_star(star_name_vector(...
            ~matches(star_name_vector, star_names)));
    end

    clear star_name_vector star_data_matrix star_importer;

    % Import satellites into simulation
    sat_importer = Sat_Data_Importer([sat_data_folder '/' sat_data_file]);
    [sat_name_vector, sat_data_matrix] = sat_importer.parse_file();
    sat_collect = Satellite_Collection(sat_name_vector, sat_data_matrix,...
        start_time);

    % If a list of satellites has been speicified in the run options,
    % remove those satellites not on the list.
    if ~isempty(sat_names)
        sat_collect.remove_sat(sat_name_vector(...
            ~matches(sat_name_vector, sat_names)));
    end
    
    clear sat_name_vector sat_data_matrix sat_importer;
    
    % Set up STK
    stk = STK_Driver(stk_version, ...
        [stk_scenario_folder '\' stk_scenario_file], ephem_step_size, ...
        occ_step_size, start_time, stop_time, oa_start_time, oa_stop_time);
    stk.start_stk(stk_visible);
    
    if import_stars
        stk.import_stars(star_collect.star_list);
    end

    if import_sats
        stk.import_sats(sat_collect.satellite_list, propagator, ...
            atmos_model);
    end
    
    toc

    % Propagate satellites in the scenario
    tic;
    stk.propagate(sat_collect.satellite_list);
    toc

    % Update the satellite ephemeris data using the propagated orbit from
    % STK
    for i = 1 : sat_collect.count()
        up_sat = sat_collect.get_sat(i);

        [time, r, v] = stk.sat_ephemeris(up_sat.name{1}, ...
            stk.analysis_start, stk.analysis_stop);

        up_sat.time_stamp = time;
        up_sat.set_icrf_position(r(:, 1), r(:, 2), r(:, 3));
        up_sat.set_icrf_velocity(v(:, 1), v(:, 2), v(:, 3));
    end

    clear i up_sat time r v;

    % Calculate occultations for each satellite and each star
    occult_calc = Occulter(start_time);

    for i = 1 : sat_collect.count()
        disp(['Satellite Number: ' num2str(i)]);
        for j = 1 : star_collect.count()
            tic;
            disp(['Star Number: ' num2str(j)]);
            a_sat = sat_collect.get_sat(i);
            a_star = star_collect.get_star(j);

            [time, r_sat, r_star] = stk.sat_star_access(a_star.name{1}, ...
                a_sat.name{1}, stk.analysis_start, stk.analysis_stop, ...
                max_star_angle, min_sun_angle, min_moon_angle);
            if ~isempty(time)
                occultations = occult_calc.occult_1(a_star.name{1}, ...
                    time, r_sat, r_star, occ_step_size);
                if ~isempty(occultations)
                    if truncate_data
                        occultations = occult_calc.truncate(...
                            occultations, lat_min, lat_max, ...
                            long_min, long_max, ...
                            alt_min, alt_max, ...
                            dt_helper.time_2_epsec(start_time), ...
                            dt_helper.time_2_epsec(stop_time));
                    end
                    a_sat.add_occ_data(occultations);
                end
            end
            toc
        end
    end
    
    clear i j a_sat a_star time r_sat r_star occultations occult_calc;
end

%% Plot Data
if plot_data
    % Create the output object
    plotter = Output(start_time);
    
    % Establish the start and stop times for plotting
    if length(oa_start_time) > 1
        start_t = dt_helper.time_2_epsec(oa_start_time);
        stop_t = dt_helper.time_2_epsec(oa_stop_time);
    else
        start_t = dt_helper.time_2_epsec(start_time);
        stop_t = dt_helper.time_2_epsec(stop_time);
    end

    % Iterate over satellites
    for j = 1 : sat_collect.count()
        sat = sat_collect.get_sat(j);

        % Iterate over region of interest use the alt_max variable. All
        % region of interest variables should be the same length
        for i = 1 : length(alt_max)
            plotter.region_2d_plot(sat, [lat_min(i) lat_max(i)], ...
                [long_min(i) long_max(i)], [alt_min(i) alt_max(i)], ...
                start_t, stop_t);
            plotter.star_avail(sat, star_collect,...
                [lat_min(i) lat_max(i)], [long_min(i) long_max(i)], ...
                [alt_min(i) alt_max(i)], start_t, stop_t);
        end

%         plotter.plot_classical(sat, start_t, stop_t);
        
        % Plot Occultation points for all stars on one graph
%         plotter.plot_occ_events_all(sat, [lat_min(1) ...
%         lat_max(1)], [long_min(1) long_max(1)], [alt_min(1) alt_max(1)],...
%         start_t, stop_t);
        plotter.plot_occ_events_all_2(sat, [lat_min(1) ...
        lat_max(1)], [long_min(1) long_max(1)], [alt_min(1) alt_max(1)],...
        start_t, stop_t);

        % Iterate over stars
        for k = 1 : star_collect.count()
            star_name = star_collect.get_star(k).name;

            try
                plotter.plot_occ_events_alt_roc(sat, ...
                    star_name, [lat_min(1) lat_max(1)], ...
                    [long_min(1) long_max(1)], [alt_min(1) alt_max(1)], ...
                    start_t, stop_t);
                plotter.plot_occ_events_duration(sat, ...
                    star_name, [lat_min(1) lat_max(1)], ...
                    [long_min(1) long_max(1)], [alt_min(1) alt_max(1)], ...
                    start_t, stop_t);
%                 plotter.plot_occ_sat_range(sat, ...
%                     star_name, [lat_min(1) lat_max(1)], ...
%                     [long_min(1) long_max(1)], [alt_min(1) alt_max(1)], ...
%                     start_t, stop_t);
                plotter.plot_occ_lst(sat, ...
                    star_name, [lat_min(1) lat_max(1)], ...
                    [long_min(1) long_max(1)], [alt_min(1) alt_max(1)], ...
                    start_t, stop_t);
            catch error_message
                disp(error_message);
            end
        end
        
        % Produce and display reports for each satellite
        report = plotter.occ_report(sat, star_collect, ...
            [lat_min(1) lat_max(1)], [long_min(1) long_max(1)], ...
            [alt_min(1) alt_max(1)], start_t, stop_t);
        disp(report);
    end

%     plotter.sat_occ_geo_plot(sat_collect.get_sat(1).occ_data(1));
%     plotter.region_2d_animate(sat_collect.get_sat(1), ...
%         star_collect.get_star(2), [lat_min(1) lat_max(1)], ...
%         [long_min(1) long_max(1)], 1, 'Mid', 3);
end

%% Save Data
if save_data
    % Generate time stamp for file saving
    dt_save_str = datestr(datetime('now'), 'yyyymmddHHMMSS');
    
    % Save important variables to a .mat file
    file_name = strcat(save_folder, '/', dt_save_str, '_', run_name);
    save(file_name, 'sat_collect', 'star_collect', 'dt_helper', '-v7.3');
end