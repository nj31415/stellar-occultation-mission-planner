%% Nicholas Jones - njones31@vt.edu
% The Output class is responsible for creating formatted graphs, visuals,
% etc to represent the results of the simulation.

classdef Output < handle
    properties (Constant)
        % https://www.mathworks.com/help/matlab/creating_plots/specify-line
        % -and-marker-appearance-in-plots.html
        marker_style = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'p', ...
            'h', '_', '|', 'v', '>', '<'};
        marker_size = 10;
    end

    properties (Access = private)
        dt_helper   % Datetime_Helper object
    end
    
    methods
        %% Constructor for an Output object. Used to access methods for
        % plotting data.
        % Inputs:
        % scenario_start    : String, formatted string representing the
        %                     scenario start time, follows format specified
        %                     in Datetime_Helper constructor.
        function obj = Output(scenario_start)
            obj.dt_helper = Datetime_Helper(scenario_start);
        end
        
        %% Function to plot the occultations for a satellite and stars 
        % within a region of interest, showing the whole globe.
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function geo_2d_plot(obj, sat, lat_lim, long_lim, alt_lim, ...
                start_time, stop_time)
            figure();
            worldmap("World");
            load coastlines;
            geoshow(coastlat, coastlon, 'Color', 'k');
            
            ax = [];
            star_legend = {};
            star_name = unique({sat.occ_data.star_name});

            for i = 1 : length(star_name)
                occultations = sat.search_occ_data_star(star_name(i));
                lla_data = double.empty(0, 4);

                % Reconcatenate data from various occultations into one
                % matrix
                for j = 1 : length(occultations)
                    idx = occultations(j).search_region_time(lat_lim, ...
                        long_lim, alt_lim, [start_time, stop_time]);
                    n = size(idx(idx == 1), 1);
                    lla_data(end + 1 : end + n, :) = ...
                        occultations(j).get_lla_data(idx);
                end

                if ~isempty(lla_data)
                    ax(end + 1) = scatterm(lla_data(:, 2), ...
                        lla_data(:, 3), obj.marker_size, lla_data(:, 4),...
                        obj.get_style(length(ax) + 1));
                    star_legend(end + 1) = star_name(i);
                    hold on;
                end
            end

            
            colormap default;
            c = colorbar;
            c.Label.String = 'Altitude, km';

            title(['Stellar Occultation Locations for ' ...
                obj.name_format(sat.name{1})]);
            legend(ax, obj.name_format(star_legend), 'Location', ...
                'northwestoutside');
        end

        %% Function to plot the occultations for a satellite and stars 
        % within a region of interest.
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function region_2d_plot(obj, sat, lat_lim, long_lim, alt_lim, ...
                start_time, stop_time)
            f = figure();
            worldmap(lat_lim, long_lim);
            load coastlines;
            geoshow(coastlat, coastlon, 'Color', 'k');
            
            ax = [];
            star_legend = {};
            star_name = unique({sat.occ_data.star_name});

            for i = 1 : length(star_name)
                occultations = sat.search_occ_data_star(star_name(i));
                lla_data = double.empty(0, 4);

                % Reconcatenate data from various occultations into one
                % matrix
                for j = 1 : length(occultations)
                    idx = occultations(j).search_region_time(lat_lim, ...
                        long_lim, alt_lim, [start_time, stop_time]);
                    n = size(idx(idx == 1), 1);
                    lla_data(end + 1 : end + n, :) = ...
                        occultations(j).get_lla_data(idx);
                end

                if ~isempty(lla_data)
                    ax(end + 1) = scatterm(lla_data(:, 2), ...
                        lla_data(:, 3), obj.marker_size, lla_data(:, 4),...
                        obj.get_style(length(ax) + 1));
                    star_legend(end + 1) = star_name(i);
                    hold on;
                end
            end

            caption = sprintf(['Stellar Occultation Locations for ' ...
                obj.name_format(sat.name{1}) '\n']);
            title(caption);
            l = legend(ax, obj.name_format(star_legend), 'Location',...
                'southoutside', 'orientation', 'horizontal');

            colormap default;
            c = colorbar('Location', 'eastoutside');
            c.Label.String = 'Altitude, km';
            
            f.Position = [100 100 750 650];
            c.Position = [0.9 0.1654 0.0328 0.7045];
            l_coord = l.Position;
            l.Position = [l_coord(1) l_coord(2) - 0.03 - l_coord(4) ...
                l_coord(3) l_coord(4)];
        end

        %% Function to animate the occultations for a satellite and stars 
        % within a region of interest. Intended for use with truncated data
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % star          : Star, star to animate occultations for
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % num_occ       : Int, number of occultations to animate
        % loc           : String, where to animate in the occultation time
        %                 series. 'Start': start animating at the first
        %                 occultation and continue for num_occ
        %                 occultations. 'Mid': animate at the
        %                 middlemost occultation with num_occ occultations 
        %                 on either side. 'End': animate the last
        %                 occultation and num_occ occultations previous.
        % lat_pad       : Float, padding to add to the map for latitude, in
        %                 degrees.
        function region_2d_animate(obj, sat, star, lat_lim, long_lim, ...
                num_occ, loc, lat_pad)
            f = figure();
            
            % Add latitude padding
            lat_lim = [lat_lim(1) - lat_pad lat_lim(2) + lat_pad];
            if lat_lim(1) < -90
                lat_lim(1) = -90;
            end
            if lat_lim(2) > 90
                lat_lim(2) = 90;
            end

            worldmap(lat_lim, long_lim);
            load coastlines;
            geoshow(coastlat, coastlon, 'Color', 'k');
            
            occultations = sat.search_occ_data_star(star.name{1});

            switch loc
                case 'Start'
                    occultations = occultations(1 : num_occ);
                case 'Mid'
                    mid_idx = round(length(occultations) / 2);
                    occultations = occultations(mid_idx - num_occ : ...
                        mid_idx + num_occ);
                case 'End'
                    occultations = occultations(end - num_occ : end);
                otherwise
                    error('Unknown location specification');
            end

            for i = 1 : length(occultations)
                occ_lla_data = occultations(i).get_lla_data([]);
                sat_icrf_data = occultations(i).get_sat_icrf_data([]);
                sat_lla_data = eci2lla(...
                    sat_icrf_data(:, 2 : 4) .* (10^3), ...
                    obj.dt_helper.epsec_2_date_vec(sat_icrf_data(:, 1)));

                ax_occ = [];
                ax_sat = [];

                for j = 1 : size(occ_lla_data, 1)
                    ax_occ(end + 1) = scatterm(occ_lla_data(j, 2), ...
                        occ_lla_data(j, 3), obj.marker_size, ...
                        occ_lla_data(j, 4), '*');
                    hold on;
                    ax_sat(end + 1) = scatterm(sat_lla_data(j, 1), ...
                        sat_lla_data(j, 2), obj.marker_size, ...
                        sat_lla_data(j, 3) .* (10^-3), 'o');

                    colormap turbo;
                    c = colorbar('Location', 'eastoutside');
                    c.Label.String = 'Altitude, km';
                    % Set colorbar limits
                    caxis([70 500]);

                    caption = sprintf(...
                        ['Time: ' ...
                        obj.dt_helper.epsec_2_fstr(occ_lla_data(j, 1)) ...
                        '\n']);
                    title(caption);

                    f.Position = [100 100 650 550];
                    c.Position = [0.9 0.1654 0.0328 0.7045];
                    
                    % Export plot to GIF
                    exportgraphics(gca, 'region_2d_animate.gif', ...
                        'Append', true);
                end
                delete(ax_occ);
                delete(ax_sat);
            end
        end

        %% Function to plot star occultation availability over time, 
        % as well as the star location in right ascension and declination.
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % star_collect  : Star Collection, used for extracting star
        %                 position data (RA, Dec)
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function star_avail(obj, sat, star_collect, lat_lim, long_lim, ...
                alt_lim, start_time, stop_time)
            figure();
            tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
            nexttile;
            star_legend = {};
            ra = double.empty();
            dec = double.empty();
            star_name = unique({sat.occ_data.star_name});

            for i = 1 : length(star_name)
                occultations = sat.search_occ_data_star(star_name(i));
                time_data = double.empty(0, 1);

                % Reconcatenate data from various occultations into one
                % matrix
                for j = 1 : length(occultations)
                    idx = occultations(j).search_region(lat_lim, ...
                        long_lim, alt_lim);
                    n = size(idx(idx == 1), 1);
                    time_data(end + 1 : end + n, :) = ...
                        occultations(j).get_time_data(idx);
                end

                if ~isempty(time_data)
                    j = length(star_legend) + 1;
                    time_data = obj.dt_helper.epsec_2_time(time_data);
                    plot(time_data, j * ones(size(time_data)), ...
                        obj.get_style(j))
                    star_legend(end + 1) = star_name(i);
                    star = star_collect.search_names(star_name(i));
                    ra(end + 1) = star.ra;
                    dec(end + 1) = star.dec;
                    hold on;
                end
            end
            
            if ~isempty(star_legend)
                title(['Star Availability Over Analysis Period for ', ...
                    obj.name_format(sat.name{1})]);
                xlabel('Time');
                ylabel('Star');
                yticks(1 : length(star_legend));
                yticklabels(obj.name_format(star_legend));
                xlim([obj.dt_helper.epsec_2_time(start_time) ...
                    obj.dt_helper.epsec_2_time(stop_time)]);
            end

            nexttile;
            for i = 1 : length(ra)
                plot(ra(i), dec(i), obj.get_style(i));
                hold on
            end
            
            title('Star Position');
            xlabel('Right Ascension, degrees');
            ylabel('Declination, degrees');
            legend(obj.name_format(star_legend), ...
                'Location', 'bestoutside');
        end

        %% Function to plot occultation events for a given satellite, star,
        % and region of interest. Further selection of the occultation
        % event can be done with the occultation index. If there are fewer
        % occultations than the index provided, that index will not be
        % plotted.
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % star_name     : String, name of the star to plot occultations for
        % idx           : Float vector, indexes of the occultation events
        %                 to plot
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        function plot_occ_events_idx(obj, sat, star_name, idx, lat_lim,...
                long_lim, alt_lim)
            figure();
            occultations = sat.search_occ_data_star(star_name);
            for i = 1 : length(idx)
                if idx(i) <= length(occultations)
                    occ_idx = occultations(idx(i))...
                        .search_region(lat_lim, long_lim, alt_lim);
                    if ~isempty(occ_idx(occ_idx == 1))
                        lla_data = occultations(idx(i))...
                            .get_lla_data(occ_idx);
                        scatter3(lla_data(:, 2), lla_data(:, 3), ...
                            lla_data(:, 4), obj.get_style(i));
                        hold on;
                    end
                end
            end

            title(['Occultations for ' obj.name_format(sat.name{1}) ...
                ' and ' obj.name_format(star_name)]);
            xlabel('Latitude, degrees');
            ylabel('Longitude degrees');
            zlabel('Altitude, km');
        end

        %% Function to plot occultation events for a given satellite, star,
        % and region of interest for a given time period
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % star_name     : String, name of the star to plot occultations for
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function plot_occ_events(obj, sat, star_name, lat_lim, long_lim,...
                alt_lim, start_time, stop_time)
            figure();

            ax = [];
            star_legend = {};

            for i = 1 : length(star_name)
                occultations = sat.search_occ_data_star(star_name{i});
                lla_data = double.empty(0, 4);

                % Reconcatenate data from various occultations into one
                % matrix
                for j = 1 : length(occultations)
                    idx = occultations(j).search_region_time(lat_lim, ...
                        long_lim, alt_lim, [start_time, stop_time]);
                    n = size(idx(idx == 1), 1);
                    lla_data(end + 1 : end + n, :) = ...
                        occultations(j).get_lla_data(idx);
                end

                if ~isempty(lla_data)
                    ax(end + 1) = scatter3(...
                        obj.dt_helper.epsec_2_time(lla_data(:, 1)), ...
                        lla_data(:, 2), lla_data(:, 3), 36, ...
                        lla_data(:, 4), obj.get_style(i));
                    star_legend(end + 1) = star_name(i);
                    hold on;
                end
            end

            c = colorbar;
            c.Label.String = 'Altitude, km';

            title(['Occultations for ' obj.name_format(sat.name{1})]);
            xlabel('Time')
            ylabel('Latitude, degrees');
            zlabel('Longitude, degrees');
            legend(ax, obj.name_format(star_legend), 'Location',...
                'northwestoutside');
            if ~isempty(star_legend)
                xlim([obj.dt_helper.epsec_2_time(start_time) ...
                    obj.dt_helper.epsec_2_time(stop_time)]);
            end
            ylim(lat_lim);
            zlim(long_lim);
        end

        %% Function to plot all occultation events for a given satellite
        % and region of interest for a given time period
        % plotted.
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function plot_occ_events_all(obj, sat, lat_lim, long_lim,...
                alt_lim, start_time, stop_time)
            figure();

            ax = [];
            star_legend = {};
            star_name = unique({sat.occ_data.star_name});

            for i = 1 : length(star_name)
                occultations = sat.search_occ_data_star(star_name(i));
                lla_data = double.empty(0, 4);

                % Reconcatenate data from various occultations into one
                % matrix
                for j = 1 : length(occultations)
                    idx = occultations(j).search_region_time(lat_lim, ...
                        long_lim, alt_lim, [start_time, stop_time]);
                    n = size(idx(idx == 1), 1);
                    lla_data(end + 1 : end + n, :) = ...
                        occultations(j).get_lla_data(idx);
                end

                if ~isempty(lla_data)
                    ax(end + 1) = scatter3(...
                        obj.dt_helper.epsec_2_time(lla_data(:, 1)), ...
                        lla_data(:, 2), lla_data(:, 3), 36, ...
                        lla_data(:, 4), obj.get_style(i));
                    star_legend(end + 1) = star_name(i);
                    hold on;
                end
            end

            c = colorbar;
            c.Label.String = 'Altitude, km';

            title(['Occultations for ' obj.name_format(sat.name{1})]);
            xlabel('Time')
            ylabel('Latitude, degrees');
            zlabel('Longitude, degrees');
            legend(ax, obj.name_format(star_legend), 'Location',...
                'northwestoutside');
            if ~isempty(star_legend)
                xlim([obj.dt_helper.epsec_2_time(start_time) ...
                    obj.dt_helper.epsec_2_time(stop_time)]);
            end
            ylim(lat_lim);
            zlim(long_lim);
        end

        %% Function to plot all occultation events for a given satellite
        % and region of interest for a given time period
        % plotted.
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function plot_occ_events_all_2(obj, sat, lat_lim, long_lim,...
                alt_lim, start_time, stop_time)
            figure();

            ax = [];
            star_legend = {};
            star_name = unique({sat.occ_data.star_name});

            for i = 1 : length(star_name)
                occultations = sat.search_occ_data_star(star_name(i));
                lla_data = double.empty(0, 4);

                % Reconcatenate data from various occultations into one
                % matrix
                for j = 1 : length(occultations)
                    idx = occultations(j).search_region_time(lat_lim, ...
                        long_lim, alt_lim, [start_time, stop_time]);
                    n = size(idx(idx == 1), 1);
                    lla_data(end + 1 : end + n, :) = ...
                        occultations(j).get_lla_data(idx);
                end

                if ~isempty(lla_data)
                    ax(end + 1) = scatter(...
                        obj.dt_helper.epsec_2_time(lla_data(:, 1)), ...
                        lla_data(:, 2), 36, lla_data(:, 4), ...
                        obj.get_style(i));
                    star_legend(end + 1) = star_name(i);
                    hold on;
                end
            end

            c = colorbar;
            c.Label.String = 'Altitude, km';

            title(['Occultations for ' obj.name_format(sat.name{1})]);
            xlabel('Time')
            ylabel('Latitude, degrees');
            legend(ax, obj.name_format(star_legend), 'Location',...
                'northwestoutside');
            if ~isempty(star_legend)
                xlim([obj.dt_helper.epsec_2_time(start_time) ...
                    obj.dt_helper.epsec_2_time(stop_time)]);
            end
            ylim(lat_lim);
        end
        
        %% Function to plot the local solar time for the measurement
        % tangent point for a given satellite, star, and region of interest
        % in a specified time
        % period
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % star_name     : String in cell array, name of the star to plot
        %                 occultations for
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function plot_occ_lst(obj, sat, star_name, lat_lim, long_lim, ...
                alt_lim, start_time, stop_time)
            occultations = sat.search_occ_data_star(star_name);
            lst_data = double.empty(0, 2);

            % Reconcatenate data from various occultations into one
            % matrix
            for j = 1 : length(occultations)
                idx = occultations(j).search_region_time(lat_lim, ...
                    long_lim, alt_lim, [start_time, stop_time]);
                n = size(idx(idx == 1), 1);
                if n > 0
                    time = occultations(j).get_time_data(idx);
                    lst = occultations(j).local_solar_time(idx);
                    lst_data(end + 1 : end + n, :) = ...
                        [time lst];
                end
            end

            if ~isempty(lst_data)
                figure();
                plot(obj.dt_helper.epsec_2_time(lst_data(:, 1)), ...
                    lst_data(:, 2), 'k*');
                hold on;

                title(['Local Solar Time of the Measurement ' ...
                    'Tangent Point for ' obj.name_format(sat.name{1}) ...
                    ' and ' obj.name_format(star_name{1})]);
                xlabel('Time');
                ylabel('Mean Local Solar Time, hr');
    
                st_time = obj.dt_helper.epsec_2_time(start_time);
                sp_time = obj.dt_helper.epsec_2_time(stop_time);
                xlim([st_time sp_time]);
                ylim([0 24]);
            end
        end

        %% Function to plot altitude rate of change for occultation events
        % for a given satellite, star, and region of interest. In a
        % specified time period.
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % star_name     : String in cell array, name of the star to plot
        %                 occultations for
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function plot_occ_events_alt_roc(obj, sat, star_name, lat_lim, ...
                long_lim, alt_lim, start_time, stop_time)
            occultations = sat.search_occ_data_star(star_name);
            alt_roc_data = double.empty(0, 2);

            % Reconcatenate data from various occultations into one
            % matrix
            for j = 1 : length(occultations)
                idx = occultations(j).search_region_time(lat_lim, ...
                    long_lim, alt_lim, [start_time, stop_time]);
                n = size(idx(idx == 1), 1);
                if n > 0
                    time = occultations(j).get_time_data(idx);
                    alt_roc = occultations(j).get_alt_roc(idx);
                    alt_roc_data(end + 1 : end + n, :) = ...
                        [time alt_roc];
                end
            end

            if ~isempty(alt_roc_data)
                figure();
                plot(obj.dt_helper.epsec_2_time(alt_roc_data(:, 1)), ...
                    alt_roc_data(:, 2), 'k*');
                hold on;

                title(['Altitude Rate of Change for Occultations for '...
                    obj.name_format(sat.name{1}) ' and ' ...
                    obj.name_format(star_name{1})]);
                xlabel('Time');
                ylabel('Altitude Rate of Change, km s^{-1}');
    
                st_time = obj.dt_helper.epsec_2_time(start_time);
                sp_time = obj.dt_helper.epsec_2_time(stop_time);
                xlim([st_time sp_time]);
            end
        end

        %% Function to plot duration of occultation events for a given
        % satellite, star, and region of interest in a specified time
        % period
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % star_name     : String in cell array, name of the star to plot
        %                 occultations for
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function plot_occ_events_duration(obj, sat, star_name, lat_lim, ...
                long_lim, alt_lim, start_time, stop_time)
            occultations = sat.search_occ_data_star(star_name);
            duration = double.empty(0, 2);

            % Reconcatenate data from various occultations into one
            % matrix
            for j = 1 : length(occultations)
                idx = occultations(j).search_region_time(lat_lim, ...
                    long_lim, alt_lim, [start_time, stop_time]);
                n = size(idx(idx == 1), 1);
                time = occultations(j).get_time_data(idx);
                if n > 0
                    duration_data = time(end) - time(1);
                    duration(end + 1, :) = [time(1) duration_data];
                end

            end

            if ~isempty(duration)
                figure();
                plot(obj.dt_helper.epsec_2_time(duration(:, 1)), ...
                    duration(:, 2), 'k*');
                hold on;

                title(['Duration of Occultations for '...
                    obj.name_format(sat.name{1}) ' and ' ...
                    obj.name_format(star_name{1})]);
                xlabel('Time');
                ylabel('Occultation Duration, s');
    
                st_time = obj.dt_helper.epsec_2_time(start_time);
                sp_time = obj.dt_helper.epsec_2_time(stop_time);
                xlim([st_time sp_time]);
            end
        end

        %% Function to plot the range between the satellite and measurement
        % tangent point for a given satellite, star, and region of interest
        % in a specified time
        % period
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % star_name     : String in cell array, name of the star to plot
        %                 occultations for
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function plot_occ_sat_range(obj, sat, star_name, lat_lim, ...
                long_lim, alt_lim, start_time, stop_time)
            occultations = sat.search_occ_data_star(star_name);
            range_data = double.empty(0, 2);

            % Reconcatenate data from various occultations into one
            % matrix
            for j = 1 : length(occultations)
                idx = occultations(j).search_region_time(lat_lim, ...
                    long_lim, alt_lim, [start_time, stop_time]);
                n = size(idx(idx == 1), 1);
                if n > 0
                    time = occultations(j).get_time_data(idx);
                    range = occultations(j).range_sat_tan(idx);
                    range_data(end + 1 : end + n, :) = ...
                        [time range];
                end
            end

            if ~isempty(range_data)
                figure();
                plot(obj.dt_helper.epsec_2_time(range_data(:, 1)), ...
                    range_data(:, 2), 'k*');
                hold on;

                title(['Range between Satellite and Measurement Tangent'...
                    ' Point for Occultations for '...
                    obj.name_format(sat.name{1}) ' and ' ...
                    obj.name_format(star_name{1})]);
                xlabel('Time');
                ylabel('Range, km');
    
                st_time = obj.dt_helper.epsec_2_time(start_time);
                sp_time = obj.dt_helper.epsec_2_time(stop_time);
                xlim([st_time sp_time]);
            end
        end

        %% Function to plot the classical orbital elements of a satellite
        % over time, extracting data from STK.
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % stk           : STK_Driver object, used to extract data.
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        % opt           : Boolean, set RAAN output option (0 => RAAN, 1 =>
        %                 LAN)
        % time_scale    : String, specifies time scale to use for the
        %                 x-axis [sec, hr, day, yr]
        function plot_classical_stk(obj, sat, stk, start_time, ...
                stop_time, opt, time_scale)
            [time, a, e, inc, RAAN, omega, nu] = ...
                stk.sat_classical(sat.name{1}, start_time, stop_time, opt);

            if strcmp(time_scale, 'hr')
                time = time ./ 3600;
            elseif strcmp(time_scale, 'day')
                time = time ./ (3600 * 24);
            elseif strcmp(time_scale, 'yr')
                time = time ./ (3600 * 24 * 365);
            end
            
            figure();
            t = tiledlayout(6, 1, 'TileSpacing', 'tight', 'Padding',...
                'tight');

            nexttile;
            plot(time, a, 'k.', 'MarkerSize', 1);
            ylabel('Semi-major Axis, km');

            nexttile;
            plot(time, e, 'k.', 'MarkerSize', 1);
            ylabel('Eccentricity');

            nexttile;
            plot(time, inc, 'k.', 'MarkerSize', 1);
            ylabel('Inclination, degrees');

            nexttile;
            plot(time, RAAN, 'k.', 'MarkerSize', 1);
            if opt == 0
                ylabel('Right Ascension of Ascending Node, degrees');
            else
                ylabel('Longitude of Ascending Node, degrees');
            end

            nexttile;
            plot(time, omega, 'k.', 'MarkerSize', 1);
            ylabel('Argument of Perigee, degrees');

            nexttile;
            plot(time, nu, 'k.', 'MarkerSize', 1);
            ylabel('True Anomaly, degrees');
            xlabel('Time, days');

            title(t, ['Classical Orbital Elements for ' ...
                obj.name_format(sat.name{1})]);
        end

        %% Function to plot the altitude of a satellite over time.
        % Inputs:
        % sat           : Satellite, satellite to plot data for.
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        function plot_classical(obj, sat, start_time, stop_time)
            ind = sat.search_time([start_time, stop_time]);
            [a, e, inc, RAAN, omega, nu] = sat.ephem_2_classical_elem(ind);
            time = obj.dt_helper.epsec_2_time(sat.time_stamp(ind));
            st_time = obj.dt_helper.epsec_2_time(start_time);
            sp_time = obj.dt_helper.epsec_2_time(stop_time);
            
            figure();
            t = tiledlayout(6, 1, 'TileSpacing', 'tight', 'Padding',...
                'tight');

            nexttile;
            plot(time, a, 'k.', 'MarkerSize', 1);
            ylabel('Semi-major Axis, km');
            xlim([st_time sp_time]);

            nexttile;
            plot(time, e, 'k.', 'MarkerSize', 1);
            ylabel('Eccentricity');
            xlim([st_time sp_time]);

            nexttile;
            plot(time, inc, 'k.', 'MarkerSize', 1);
            ylabel('Inclination, degrees');
            xlim([st_time sp_time]);

            nexttile;
            plot(time, RAAN, 'k.', 'MarkerSize', 1);
            ylabel('Right Ascension of Ascending Node, degrees');
            xlim([st_time sp_time]);

            nexttile;
            plot(time, omega, 'k.', 'MarkerSize', 1);
            ylabel('Argument of Perigee, degrees');
            xlim([st_time sp_time]);

            nexttile;
            plot(time, nu, 'k.', 'MarkerSize', 1);
            ylabel('True Anomaly, degrees');
            xlabel('Time (Days)');
            xlim([st_time sp_time]);

            title(t, ['Classical Orbital Elements for ' ...
                obj.name_format(sat.name{1})]);
        end

        %% Function to plot the latitude, longitude, and altitude of the
        % tangent point and satellite for an occultation.
        % Inputs:
        % occ_event     : Occ_Event, the occultation to plot
        function sat_occ_geo_plot(obj, occ_event)
            occ_lla_data = occ_event.get_lla_data([]);
            sat_icrf_data = occ_event.get_sat_icrf_data([]);
            sat_lla_data = eci2lla(sat_icrf_data(:, 2 : 4) .* (10^3), ...
                obj.dt_helper.epsec_2_date_vec(sat_icrf_data(:, 1)));

            figure();
            g = geoglobe(uifigure);
            g.Terrain = 'none';
            p1 = geoplot3(g, occ_lla_data(:, 2), occ_lla_data(:, 3), ...
                occ_lla_data(:, 4) .* (10^3), 'mo');
            p1.HeightReference = 'ellipsoid';
            hold(g, 'on');
            p2 = geoplot3(g, sat_lla_data(:, 1), sat_lla_data(:, 2), ...
                sat_lla_data(:, 3), 'ro');
            p2.HeightReference = 'ellipsoid';
        end

        %% Function to produce an occultation report for a given satellite.
        % Inputs:
        % sat           : Satellite, the satellite to produce the report
        %                 for
        % star_collect  : Star Collection, contains star data for the
        %                 simulation
        % lat_lim       : Float vector, specification of the min and max
        %                 latitude of the region of interest in deg.
        %                 [min_lat max_lat]
        % long_lim      : Float vector, specification of the min and max
        %                 longitude of the region of interest in deg. 
        %                 [min_long max_long]
        % alt_lim       : Float vector, specification of the min and max
        %                 altitude of the region of interest in km.
        %                 [min_alt max_alt]
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        % Outputs:
        % report_array  : Cell array, summary table for occultations
        function report_array = occ_report(obj, sat, star_collect, ...
                lat_lim, long_lim, alt_lim, start_time, stop_time)
            star_name = unique({sat.occ_data.star_name});

            if isempty(star_name)
                report_array = 'No Occultations';
            else

                star_list = Star.empty();
                
                num_occ = double.empty();
                max_time = NaT('TimeZone', 'UTC');
                min_time = NaT('TimeZone', 'UTC');
                max_lat = double.empty();
                min_lat = double.empty();
                max_alt = double.empty();
                min_alt = double.empty();
                max_alt_roc = double.empty();
    
                val_idx = 1;
    
                for i = 1 : length(star_name)
                    occultations = sat.search_occ_data_star(star_name(i));
                    lla_data = double.empty(0, 4);
                    alt_roc_data = double.empty();
                    occ_count = 0;
    
                    % Reconcatenate data from various occultations into one
                    % matrix
                    for j = 1 : length(occultations)
                        idx = occultations(j).search_region_time(...
                            lat_lim, long_lim, alt_lim, [start_time, ...
                            stop_time]);
                        n = size(idx(idx == 1), 1);
                        lla_data(end + 1 : end + n, :) = ...
                            occultations(j).get_lla_data(idx);
    
                        if n > 0
                            alt_roc = occultations(j).get_alt_roc(idx);
                            alt_roc_data(end + 1 : end + n, :) = alt_roc;
                            occ_count = occ_count + 1;
                        end
                    end
    
                    if ~isempty(lla_data)
                        star_list(val_idx) = star_collect.search_names(...
                            star_name(i));
                        num_occ(val_idx) = occ_count;
                        max_time(val_idx) = obj.dt_helper.epsec_2_time(...
                            max(lla_data(:, 1)));
                        min_time(val_idx) = obj.dt_helper.epsec_2_time(...
                            min(lla_data(:, 1)));
                        max_lat(val_idx) = max(lla_data(:, 2));
                        min_lat(val_idx) = min(lla_data(:, 2));
                        max_alt(val_idx) = max(lla_data(:, 4));
                        min_alt(val_idx) = min(lla_data(:, 4));
                        max_alt_roc(val_idx) = max(abs(alt_roc_data));
    
                        val_idx = val_idx + 1;
                    end
                end
                
                star_name = convertCharsToStrings([star_list.name])';
                hd_num = [star_list.hd_num]';
                ra = [star_list.ra]';           % deg - right ascension
                dec = [star_list.dec]';         % deg - declination
                num_occ = num_occ';             % number of occultations
                min_time = min_time';           % First occultation time
                max_time = max_time';           % Last occultation time
                avail_dur = max_time - min_time;
                avail_dur.Format = 'd';         % days - Availability
                                                % duration
                min_lat = min_lat';             % deg - latitude
                max_lat = max_lat';             % deg - latitude
                delta_lat = max_lat - min_lat;  % deg - delta latitude
                min_alt = min_alt';             % km - altitude
                max_alt = max_alt';             % km - altitude
                max_alt_roc = max_alt_roc';     % km s^-1 - altitude roc
    
                report_array = table(star_name, hd_num, ra, dec, ...
                    num_occ, min_time, max_time, avail_dur, min_lat, ...
                    max_lat, delta_lat, min_alt, max_alt, max_alt_roc);
            end
        end

        %% Function to get a marker style based on the supplied index.
        % Inputs:
        % idx   : Int, style number, will be modded to access valid index
        %         in marker_style vector
        % Outputs:
        % m_s   : character vector, the marker style to use
        function m_s = get_style(obj, idx)
            m_s = obj.marker_style(mod(idx, length(obj.marker_style)) + 1);
            m_s = m_s{:};
        end

        %% Function to format satellite and star names to display properly
        % in plots and legends
        % Inputs:
        % str   : String or Cell Array, strings to remove '_' and replace
        %         with 'a'
        function form_str = name_format(~, str)
            form_str = strrep(str, '_', ' ');
        end
    end
end

