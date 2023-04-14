%% Nicholas Jones - njones31@vt.edu
% STK_Driver
% The STK_Driver class controls interfacing between the MATLAB script and
% STK for creating satellite and star objects in the STK scenario,
% extracting data from those objects, calculating access windows, and
% propagating along the start and stop times of the scenario.

classdef STK_Driver < handle
    properties (SetAccess = immutable)
        stk_version = 12;       % Default value
        scenario_file
        ephem_time_step = 10;   % Default value - sec
        occ_time_step = 5;      % Default value - sec    
        prop_lim = 100;         % km - propagation altitude limit for HPOP
    end

    properties (Access = private)
        dt_helper           % Datetime_Helper object
        uiApplication       % STK application handle
        root                % STK application root
        scenario            % STK scenario object
    end

    properties (SetAccess = public)
        scenario_start
        scenario_stop
        analysis_start
        analysis_stop
    end
    
    methods
        %% Constructor for STK_Driver class. Used to populate STK_Driver
        % properties before starting STK.
        % Inputs:
        % stk_version       : int, STK Version number. Should be 11 or 12.
        %                     Note all testing performed with STK 12.
        % scenario_file     : String, file path to scenario file for use in
        %                     the analysis.
        % ephem_time_step   : float, time step for retrieving satellite
        %                     ephemeris in secs.
        % occ_time_step     : float, time step for determining 
        %                     satellite-star accesses in sec
        % scenario_start    : String, time specifiation of the scenario
        %                     start time as a string. Formatted as:
        %                     yyyy:mm:dd:HH:MM:SS.FFF'
        % scenario_stop     : String or float, time specifiation of the
        %                     scenario stop time as a string or a float in
        %                     Epsec.
        % analysis_start    : String or float, time specifiation of the
        %                     analysis start time as a string or a float in
        %                     Epsec. Optional. If unspecified, analysis
        %                     will begin at the scenario start time.
        % analysis_stop     : String or float, time specifiation of the
        %                     analysis stop time as a string or a float in
        %                     Epsec. Optional. If unspecified, analysis
        %                     will stop at the scenario stop time.
        function obj = STK_Driver(stk_version, scenario_file, ...
                ephem_time_step, occ_time_step, scenario_start, ...
                scenario_stop, analysis_start, ...
                analysis_stop)
            obj.stk_version = stk_version;

            if isfile(scenario_file)
                obj.scenario_file = scenario_file;
            else
                error(['STK_Driver::set.scenario_file: Could not find', ...
                    ' the specified file']);
            end

            obj.ephem_time_step = ephem_time_step;
            obj.occ_time_step = occ_time_step;
            obj.scenario_start = scenario_start;

            obj.dt_helper = Datetime_Helper(scenario_start);

            obj.scenario_stop = scenario_stop;

            if exist('analysis_start', 'var') && ~isempty(analysis_start)
                obj.analysis_start = analysis_start;
            else
                obj.analysis_start = scenario_start;
            end

            if exist('analysis_stop', 'var') && ~isempty(analysis_stop)
                obj.analysis_stop = analysis_stop;
            else
                obj.analysis_stop = scenario_stop;
            end
        end
        
        %% Get the STK Version number.
        % Outputs:
        % stk_version   : int, STK version number
        function stk_version = get.stk_version(obj)
            stk_version = obj.stk_version;
        end

        %% Get the scenario file path.
        % Outputs:
        % scenario_file : String, file path to the STK scenario
        function scenario_file = get.scenario_file(obj)
            scenario_file = obj.scenario_file;
        end

        %% Get the STK ephemeris time step in sec
        % Outputs:
        % time_step : Float, STK ephemeris time step
        function time_step = get.ephem_time_step(obj)
            time_step = obj.ephem_time_step;
        end

        %% Get the STK occultation time step in sec
        % Outputs:
        % time_step : Float, STK occultation time step
        function time_step = get.occ_time_step(obj)
            time_step = obj.occ_time_step;
        end

        %% Get the STK scenario start time as a formatted string
        % Outputs:
        % scenario_start    : String, formatted string representing the
        %                     scenario start time.
        function scenario_start = get.scenario_start(obj)
            scenario_start = obj.scenario_start;
        end

        %% Set the STK scenario start time as a formatted string
        % Inputs:
        % scenario_start    : String, formatted string representing the
        %                     scenario start time
        function set.scenario_start(obj, scenario_start)
            obj.scenario_start = scenario_start;
        end

        %% Get the STK scenario stop time in EpSec
        % Outputs:
        % scenario_stop : Float, scenario stop time in EpSec.
        function scenario_stop = get.scenario_stop(obj)
            scenario_stop = obj.scenario_stop;
        end

        %% Set the STK scenario stop time.
        % Inputs:
        % scenario_stop : String or float, scenario stop time as a
        %                 formatted time string (same format as
        %                 scenario_start string), or as a value in EpSec.
        function set.scenario_stop(obj, scenario_stop)
            if isa(scenario_stop, 'string') || isa(scenario_stop, 'char')
                obj.scenario_stop = obj.dt_helper...
                    .time_2_epsec(scenario_stop);
            elseif isa(scenario_stop, 'double')
                obj.scenario_stop = scenario_stop;
            else
                error('STK_Driver::set.scenario_stop: Unknown format.');
            end
        end
        
        %% Get the analysis start time in EpSec
        % Outputs:
        % analysis_start    : Float, analysis start time in EpSec.
        function analysis_start = get.analysis_start(obj)
            analysis_start = obj.analysis_start;
        end

        %% Set the analysis start time
        % Inputs:
        % analysis_start    : String or float, analysis start time as a
        %                     formatted time string (same format as
        %                     scenario_start string), or as a value in 
        %                     EpSec.
        function set.analysis_start(obj, analysis_start)
            if isa(analysis_start, 'string') || isa(analysis_start, 'char')
                obj.analysis_start = obj.dt_helper...
                    .time_2_epsec(analysis_start);
            elseif isa(analysis_start, 'double')
                obj.analysis_start = analysis_start;
            else
                error('STK_Driver::set.analysis_start: Unknown format.');
            end
        end

        %% Get the analysis stop time in EpSec
        % Outputs:
        % analysis_stop : Float, analysis stop time in EpSec.
        function analysis_stop = get.analysis_stop(obj)
            analysis_stop = obj.analysis_stop;
        end

        %% Set the analysis stop time
        % Inputs:
        % analysis_stop : String or float, analysis stop time as a
        %                 formatted time string (same format as
        %                 scenario_start string), or as a value in EpSec.
        function set.analysis_stop(obj, analysis_stop)
            if isa(analysis_stop, 'string') || isa(analysis_stop, 'char')
                obj.analysis_stop = obj.dt_helper...
                    .time_2_epsec(analysis_stop);
            elseif isa(analysis_stop, 'double')
                obj.analysis_stop = analysis_stop;
            else
                error('STK_Driver::set.analysis_stop: Unknown format.');
            end
        end

        %% Start STK application
        % Inputs:
        % visible   : Boolean, controls whether STK window is visible
        function start_stk(obj, visible)
            ver_string = ['STK' num2str(obj.stk_version) '.Application'];

            try
                obj.uiApplication = actxserver(ver_string);
                obj.uiApplication.UserControl = 1;
                obj.uiApplication.Visible = visible;
                obj.root = obj.uiApplication.Personality2;
                obj.root.LoadScenario([obj.scenario_file]);
                obj.scenario = obj.root.CurrentScenario;
                
                % Set the STK date format to accept absolute time in a
                % formatted string.
                obj.root.UnitPreferences.Item('DateFormat')...
                    .SetCurrentUnit('YYYY:MM:DD');

                obj.scenario.SetTimePeriod(obj.scenario_start, ...
                    obj.dt_helper.epsec_2_fstr(obj.scenario_stop));
                obj.scenario.Epoch = obj.scenario_start;
                obj.root.Rewind();
                
                % Set the STK date format to EpSec. All calculations are
                % carried out in this unit.
                obj.root.UnitPreferences.Item('DateFormat')...
                    .SetCurrentUnit('EpSec');
            catch error_message
                % Close the STK application if there is an error
                % encountered during start up.
                obj.close_stk();
                disp(error_message);
                disp(error_message.stack(1).file);
                disp(error_message.stack(1).name);
                disp(error_message.stack(1).line);
            end
        end

        %% Close STK
        function close_stk(obj)
            obj.root.CloseScenario();
            obj.uiApplication.Quit();
        end

        %% Function to import stars into the STK scenario
        % Inputs:
        % star_vector   : Star object vector, vector of star objects,
        %                 ideally from a Star_Collection
        function import_stars(obj, star_vector)
            % Unit preferences are changed to match the incoming data
            obj.root.UnitPreferences.Item('Time').SetCurrentUnit('yr');
            
            % Create new stars and assign properties
            for i = 1 : size(star_vector, 1)
                new_star = obj.scenario.Children.New('eStar', ...
                    star_vector(i).name{1});

                obj.root.UnitPreferences.Item('Angle')...
                    .SetCurrentUnit('deg');
                new_star.LocationRightAscension = star_vector(i).ra;
                new_star.LocationDeclination = star_vector(i).dec;

                obj.root.UnitPreferences.Item('Angle')...
                    .SetCurrentUnit('arcSec');
                new_star.ProperMotionRightAscension = star_vector(i)...
                    .ra_proper_motion;
                new_star.ProperMotionDeclination = star_vector(i)...
                    .dec_proper_motion;
                new_star.Parallax = star_vector(i).parallax;
            end
            
            % Rest unit preferences to default
            obj.root.UnitPreferences.Item('Time').SetCurrentUnit('sec');
            obj.root.UnitPreferences.Item('Angle').SetCurrentUnit('deg');
        end

        %% Function to import satellites into STK scenario
        % Inputs:
        % sat_vector    : Satellite object vector, vector of satellite
        %                 objects, ideally from a Satellite_Collection
        % propagator    : String, specifies the propagator to use for each
        % atmos_model   : String, specifies the atmospheric model to use
        function import_sats(obj, sat_vector, propagator, atmos_model)
            % Create new satelllites and assign properties
            for i = 1 : size(sat_vector, 1)
                new_sat = obj.scenario.Children...
                    .New('eSatellite', sat_vector(i).name{1});
                
                % Assign satellite propagator
                if strcmp(propagator, 'TwoBody')
                    new_sat.SetPropagatorType('ePropagatorTwoBody');
                elseif strcmp(propagator, 'J2')
                    new_sat.SetPropagatorType('ePropagatorJ2Perturbation');
                elseif strcmp(propagator, 'HPOP')
                    new_sat.SetPropagatorType('ePropagatorHPOP');
                    obj.set_hpop_opt(new_sat, sat_vector(i).drag_coeff, ...
                        sat_vector(i).area_mass_ratio, atmos_model);
                else
                    error(['STK_Driver: Unrecognized propagator in ', ...
                        'run options']);
                end

                % Retrieve and update initial state of satellite
                keplerian = new_sat.Propagator.InitialState...
                    .Representation.ConvertTo('eOrbitStateClassical');
                keplerian.SizeShapeType = 'eSizeShapeSemimajorAxis';

                keplerian.SizeShape.SemiMajorAxis = sat_vector(i).a_init;
                keplerian.SizeShape.Eccentricity = sat_vector(i).e_init;

                keplerian.Orientation.Inclination = sat_vector(i).i_init;

                keplerian.Orientation.AscNodeType = 'eAscNodeLAN';
                keplerian.Orientation.AscNode.Value = ...
                    sat_vector(i).lan_init;

                keplerian.Orientation.ArgOfPerigee = ...
                    sat_vector(i).omega_init;

                keplerian.LocationType = 'eLocationTrueAnomaly';
                keplerian.Location.Value = sat_vector(i).nu_init;

                new_sat.Propagator.InitialState.Representation...
                    .Assign(keplerian);
            end
        end

        %% Function to set HPOP propagator properties
        % Inputs:
        % sat               : STK satellite object, Satellite to set HPOP
        %                     properties for
        % drag_coeff        : Satellite drag coefficient
        % area_mass_ratio   : Float, area to mass ratio in m^2 kg^-1
        % atmos_model       : String, formatted string for setting the
        %                     atmospheric model.
        function set_hpop_opt(obj, sat, drag_coeff, area_mass_ratio, ...
                atmos_model)
            hpop_prop = sat.Propagator;
            hpop_force_model = hpop_prop.ForceModel;
            hpop_force_model.Drag.Use = true;
            hpop_drag = hpop_force_model.Drag.DragModel;
            hpop_drag.AreaMassRatio = area_mass_ratio;
            hpop_drag.Cd = drag_coeff;
            hpop_force_model.Drag.AtmosphericDensityModel = atmos_model;

            % Disable solar radiation pressure calculations
            hpop_force_model.SolarRadiationPressure.Use = 0;

            % Set lower limit of propagation to 100 km altitude
            hpop_prop.Integrator.DoNotPropagateBelowAlt = obj.prop_lim;
        end

        %% Function to propagate the satellites in the scenario
        function propagate(obj, sat_list)
            for i = 1 : length(sat_list)
                sat = obj.scenario.Children.Item(sat_list(i).name{1});
                sat.Propagator.Propagate();
            end
        end

        %% Function to extract star location in ICRF coordinates from the
        % STK scenario. Exported as column vectors
        % Inputs:
        % name          : String, the name of the star to extract the data
        %                 for.
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        % Outputs:
        % time          : Float vector, time stamps for the ephemeris data
        %                 in EpSec.
        % r             : Float vector, position data in ICRF frame
        function [time, r] = star_ephemeris(obj, name, start_time, ...
                stop_time)
            rpt_elems = {'Time'; 'x'; 'y'; 'z'};

            star = obj.scenario.Children.Item(name);

            star_dp = star.DataProviders...
                .GetDataPrvTimeVarFromPath('Vectors(ICRF)//Location');
            star_data = star_dp.ExecElements(start_time, stop_time, ...
                obj.ephem_time_step, rpt_elems);
            
            star_time = star_data.DataSets.Item(0).GetValues();

            sz = size(star_time, 1);
            star_pos = cell(sz, 3);

            star_pos(:, 1) = star_data.DataSets.Item(1).GetValues();
            star_pos(:, 2) = star_data.DataSets.Item(2).GetValues();
            star_pos(:, 3) = star_data.DataSets.Item(3).GetValues();

            time = cell2mat(star_time);
            r = cell2mat(star_pos);
        end

        %% Function to extract satellite position and velocity in ICRF
        % coordinates from the STK scenario. Exported as column vectors
        % Inputs:
        % name          : String, the name of the satellite to extract data
        %                 for.
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        % Outputs:
        % time          : Float vector, time stamps for the ephemeris data
        %                 in EpSec.
        % r             : Float vector, position data in ICRF frame
        % v             : Float vector, velocity data in ICRF frame
        function [time, r, v] = sat_ephemeris(obj, name, start_time, ...
                stop_time)
            rpt_elems = {'Time'; 'x'; 'y'; 'z';};

            sat = obj.scenario.Children.Item(name);

            sat_pos_dp = sat.DataProviders...
                .GetDataPrvTimeVarFromPath('Cartesian Position//ICRF');
            sat_vel_dp = sat.DataProviders...
                .GetDataPrvTimeVarFromPath('Cartesian Velocity//ICRF');

            sat_pos_data = sat_pos_dp.ExecElements(start_time,...
                stop_time, obj.ephem_time_step, rpt_elems);
            sat_vel_data = sat_vel_dp.ExecElements(start_time,...
                stop_time, obj.ephem_time_step, rpt_elems);

            sat_time = sat_pos_data.DataSets.Item(0).GetValues();

            sz = size(sat_time, 1);
            sat_pos = cell(sz, 3);
            sat_vel = cell(sz, 3);

            sat_pos(:, 1) = sat_pos_data.DataSets.Item(1).GetValues();
            sat_pos(:, 2) = sat_pos_data.DataSets.Item(2).GetValues();
            sat_pos(:, 3) = sat_pos_data.DataSets.Item(3).GetValues();

            sat_vel(:, 1) = sat_vel_data.DataSets.Item(1).GetValues();
            sat_vel(:, 2) = sat_vel_data.DataSets.Item(2).GetValues();
            sat_vel(:, 3) = sat_vel_data.DataSets.Item(3).GetValues();

            time = cell2mat(sat_time);
            r = cell2mat(sat_pos);
            v = cell2mat(sat_vel);
        end

        %% Function to get ephemeris data at access periods between a
        % target star and satellite.
        % Inputs:
        % star_name         : String, name of the target star
        % sat_name          : String, name of the targetting satellite
        % start_time        : Float, extraction start time in EpSec, must
        %                     be within the scenario start and stop times
        %                     ( >= 0, < stop_time)
        % stop_time         : Float, extraction stop time in EpSec, must be
        %                     within the scenario start and stop times 
        %                     ( > start_time, <= scenario_stop)
        % max_star_angle    : Float, maximum angle between the star and
        %                     orbit plane, specified in degrees.
        % min_sun_angle     : Float, minimum angle between the instrument
        %                     (vector from spacecraft to star) and the Sun
        %                     to be considered a valid occultation.
        % min_moon_angle    : Float, minimum angle between the instrument
        %                     (vector from spacecraft to star) and the Moon
        %                     to be considered a valid occultation.
        % Outputs:
        % time              : Float vector, time stamps for access times in
        %                     EpSec.
        % r_sat             : Float vector, sat position data in ICRF frame
        %                     in km
        % r_star            : Float vector, star position data in ICRF
        %                     frame
        function [time, r_sat, r_star] = sat_star_access(obj, star_name,...
                sat_name, start_time, stop_time, max_star_angle, ...
                min_sun_angle, min_moon_angle)
            rpt_elems_from = {'Time'; 'x'; 'y'; 'z'};
            rpt_elems_to = {'Time'; 'x/Range'; 'y/Range'; 'z/Range'};
            
            % Retrieve references to satellite and star STK objects
            star = obj.scenario.Children.Item(star_name);
            sat = obj.scenario.Children.Item(sat_name);
            
            obj.root.BeginUpdate();
            % Set up access constraints for star angle from orbit plane,
            % and Sun and Moon keep-outs
            vgt_sat = sat.Vgt;
            vgt_star = star.Vgt;

            orb_normal = vgt_sat.Vectors.Item('Orbit_Normal');
            sat_sun = vgt_sat.Vectors.Item('Sun');
            sat_moon = vgt_sat.Vectors.Item('Moon');
            star_pos = vgt_star.Vectors.Item('Location');

            angle_factory = vgt_sat.Angles.Factory;
            
            % Set up angle between orbit normal and star position to act as
            % access constraint to enforce max_star_angle
            star_orb_angle = angle_factory.Create(...
                ['star_orbit_normal_angle_' star_name], ...
                'Angle between orbit normal and star position', ...
                'eCrdnAngleTypeBetweenVectors');
            star_orb_angle.FromVector.SetVector(orb_normal);
            star_orb_angle.ToVector.SetVector(star_pos);

            % Set up angle between sun and star position to act as access
            % constraint to enforce min_sun_angle
            sun_star_angle = angle_factory.Create(...
                ['sun_star_angle_' star_name], ...
                'Angle between the sun and star position', ...
                'eCrdnAngleTypeBetweenVectors');
            sun_star_angle.FromVector.SetVector(sat_sun);
            sun_star_angle.ToVector.SetVector(star_pos);

            % Set up angle between moon and star position to act as access
            % constraint to enforce min_moon_angle
            moon_star_angle = angle_factory.Create(...
                ['moon_star_angle_' star_name], ...
                'Angle between the moon and star position', ...
                'eCrdnAngleTypeBetweenVectors');
            moon_star_angle.FromVector.SetVector(sat_moon);
            moon_star_angle.ToVector.SetVector(star_pos);
            
            % Enable max_star_angle constraint
            awb_access_constraints = sat.AccessConstraints.AWBConstraints;
            constraint = awb_access_constraints.AddConstraint(...
                'eCstrAWBAngle', star_orb_angle.QualifiedPath);
            constraint.EnableMin = true;
            constraint.Min = 90 - max_star_angle;
            constraint.EnableMax = true;
            constraint.Max = 90 + max_star_angle;

            % Enable min_sun_angle constraint
            constraint = awb_access_constraints.AddConstraint(...
                'eCstrAWBAngle', sun_star_angle.QualifiedPath);
            constraint.EnableMin = true;
            constraint.Min = min_sun_angle;

             % Enable min_moon_angle constraint
            constraint = awb_access_constraints.AddConstraint(...
                'eCstrAWBAngle', moon_star_angle.QualifiedPath);
            constraint.EnableMin = true;
            constraint.Min = min_moon_angle;

            obj.root.EndUpdate();

            % Calculate access between the satellite and the star
            access = sat.GetAccessToObject(star);

            access.AccessTimePeriod = 'eUserSpecAccessTime';
            access.SpecifyAccessTimePeriod(start_time, stop_time);

            access.ComputeAccess();
            
            % Extract ephemeris data during access periods
            access_sat_dp = access.DataProviders...
                .Item('From Position Velocity').Group().Item('ICRF');
            access_star_dp = access.DataProviders...
                .Item('To Position Velocity').Group().Item('ICRF');

            access_sat_data = access_sat_dp.ExecElements(start_time, ...
                stop_time, obj.occ_time_step, rpt_elems_from);
            access_star_data = access_star_dp.ExecElements(start_time, ...
                stop_time, obj.occ_time_step, rpt_elems_to);
            
            % Grow arrays as access periods are iterated through
            time = double.empty(0, 1);
            r_sat = double.empty(0, 3);
            r_star = double.empty(0, 3);

            for k = 1 : 4 : access_sat_data.DataSets.Count
                t_stk = cell2mat(access_sat_data.DataSets.Item(k - 1)...
                    .GetValues());
                n = size(t_stk, 1);
                time(end + 1 : end + n) = t_stk;

                r_sat(end + 1 : end + n, 1) = cell2mat(...
                    access_sat_data.DataSets.Item(k).GetValues());
                r_sat(end - n + 1 : end, 2) = cell2mat(...
                    access_sat_data.DataSets.Item(k + 1).GetValues());
                r_sat(end - n + 1 : end, 3) = cell2mat(...
                    access_sat_data.DataSets.Item(k + 2).GetValues());

                r_star(end + 1 : end + n, 1) = cell2mat(...
                    access_star_data.DataSets.Item(k).GetValues());
                r_star(end - n + 1 : end, 2) = cell2mat(...
                    access_star_data.DataSets.Item(k + 1).GetValues());
                r_star(end - n + 1 : end, 3) = cell2mat(...
                    access_star_data.DataSets.Item(k + 2).GetValues());
            end

            % Remove the access constraints
            obj.root.BeginUpdate();
            awb_access_constraints.RemoveConstraint('eCstrAWBAngle', ...
                star_orb_angle.QualifiedPath);
            awb_access_constraints.RemoveConstraint('eCstrAWBAngle', ...
                sun_star_angle.QualifiedPath);
            awb_access_constraints.RemoveConstraint('eCstrAWBAngle', ...
                moon_star_angle.QualifiedPath);
            obj.root.EndUpdate();
        end

        %% Function to retrieve LLA data for a satellite from STK
        % Inputs:
        % name          : String, the name of the satellite to extract data
        %                 for.
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        % Outputs:
        % time          : Float vector, time stamps for access times in
        %                 EpSec.
        % lat           : Float vector, satellite latitude in degrees
        % long          : Float vector, satellite longitude in degrees
        % alt           : Float vector, satellite altitude in degrees
        function [time, lat, long, alt] = sat_lla(obj, name, start_time,...
                stop_time)
            rpt_elems = {'Time'; 'Lat'; 'Lon'; 'Alt'};

            sat = obj.scenario.Children.Item(name);

            sat_lla_dp = sat.DataProviders.GetDataPrvTimeVarFromPath(...
                'LLA State//TrueOfDateRotating');
            sat_lla_data = sat_lla_dp.ExecElements(start_time, ...
                stop_time, obj.ephem_time_step, rpt_elems);

            sat_time = sat_lla_data.DataSets.Item(0).GetValues();
            
            sz = size(sat_time, 1);
            sat_lat = cell(sz, 1);
            sat_long = cell(sz, 1);
            sat_alt = cell(sz, 1);

            sat_lat(:) = sat_lla_data.DataSets.Item(1).GetValues();
            sat_long(:) = sat_lla_data.DataSets.Item(2).GetValues();
            sat_alt(:) = sat_lla_data.DataSets.Item(3).GetValues();

            time = cell2mat(sat_time);
            lat = cell2mat(sat_lat);
            long = cell2mat(sat_long);
            alt = cell2mat(sat_alt);
        end

        %% Function to retrieve satellite classical orbital elements
        % Inputs:
        % name          : String, the name of the satellite to extract data
        %                 for.
        % start_time    : Float, extraction start time in EpSec, must be
        %                 within the scenario start and stop times
        %                 ( >= 0, < stop_time)
        % stop_time     : Float, extraction stop time in EpSec, must be
        %                 within the scenario start and stop times 
        %                 ( > start_time, <= scenario_stop)
        % opt           : Boolean, set RAAN output option (0 => RAAN, 1 =>
        %                 LAN)
        % Outputs:
        % time          : Float vector, time stamps for simulation times in
        %                 EpSec.
        % a             : Float vector, semimajor axis in km
        % e             : Float vector, eccentricity
        % inc           : Float vector, inclination in degrees
        % RAAN          : Float vector, Right Ascension of Ascending Node
        %                 in deg, unless opt is set to 1, then LAN in
        %                 degrees
        % omega         : Float vector, argument of perapsis in degrees
        % nu            : Float vector, true anomaly in degrees
        function [time, a, e, inc, RAAN, omega, nu] = ...
                sat_classical(obj, name, start_time, stop_time, opt)
            if opt == 0
                rpt_elems = {'Time'; 'Semi-major Axis'; 'Eccentricity';...
                    'Inclination'; 'RAAN'; 'Arg of Perigee'; ...
                    'True Anomaly'};
            else
                rpt_elems = {'Time'; 'Semi-major Axis'; 'Eccentricity';...
                    'Inclination'; 'Lon Ascn Node'; 'Arg of Perigee'; ...
                    'True Anomaly'};
            end

            sat = obj.scenario.Children.Item(name);

            sat_classical_dp = sat.DataProviders...
                .GetDataPrvTimeVarFromPath('Classical Elements//ICRF');
            sat_clasical_data = sat_classical_dp.ExecElements(...
                start_time, stop_time, obj.ephem_time_step, rpt_elems);

            sat_time = sat_clasical_data.DataSets.Item(0).GetValues();

            sz = size(sat_time, 1);
            sat_a = cell(sz, 1);
            sat_e = cell(sz, 1);
            sat_inc = cell(sz, 1);
            sat_RAAN = cell(sz, 1);
            sat_omega = cell(sz, 1);
            sat_nu = cell(sz, 1);
            
            % STK returns data inconsistent with the order given in report
            % elements. Check the order in the Report and Graph Manager
            % tool for the actual return order
            if opt == 0
                sat_a(:) = sat_clasical_data.DataSets.Item(1).GetValues();
                sat_e(:) = sat_clasical_data.DataSets.Item(2).GetValues();
                sat_inc(:) = sat_clasical_data.DataSets.Item(3)...
                    .GetValues();
                sat_RAAN(:) = sat_clasical_data.DataSets.Item(4)...
                    .GetValues();
                sat_omega(:) = sat_clasical_data.DataSets.Item(5)...
                    .GetValues();
                sat_nu(:) = sat_clasical_data.DataSets.Item(6).GetValues();
            else
                sat_a(:) = sat_clasical_data.DataSets.Item(1).GetValues();
                sat_e(:) = sat_clasical_data.DataSets.Item(2).GetValues();
                sat_inc(:) = sat_clasical_data.DataSets.Item(3)...
                    .GetValues();
                sat_RAAN(:) = sat_clasical_data.DataSets.Item(6)...
                    .GetValues();
                sat_omega(:) = sat_clasical_data.DataSets.Item(4)...
                    .GetValues();
                sat_nu(:) = sat_clasical_data.DataSets.Item(5).GetValues();
            end

            time = cell2mat(sat_time);
            a = cell2mat(sat_a);
            e = cell2mat(sat_e);
            inc = cell2mat(sat_inc);
            RAAN = cell2mat(sat_RAAN);
            omega = cell2mat(sat_omega);
            nu = cell2mat(sat_nu);
        end
    end
end