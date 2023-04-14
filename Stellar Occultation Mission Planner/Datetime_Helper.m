%% Nicholas Jones - njones31@vt.edu
% Datetime_Helper
% The Datetime_Helper class assists as a front-end / back-end interface,
% translating user-input from date-time strings to EpSecs, as well as any
% translations that may need to occur on the back-end for different
% function inputs. All times should be in UTC where applicable.

classdef Datetime_Helper < handle
    properties (Constant)
        % Format specifier for input and output time strings
        format_in = 'yyyy:MM:dd:HH:mm:ss.SSS';
        format_out = 'yyyy:mm:dd:HH:MM:SS.FFF';
    end

    properties (SetAccess = public)
        init_time
    end

    properties (SetAccess = private)
        init_time_datetime
    end
    
    methods
        %% Constructor for Datetime_Helper object.
        % Inputs:
        % init_time : String, formatted sting representing the scenario
        %             start time, follows format_in.
        function obj = Datetime_Helper(init_time)
            obj.init_time = init_time;
        end
        
        %% Get the initial time string.
        % Outputs:
        % init_time : String, formatted sting representing the scenario
        %             start time, follows format_in.
        function init_time = get.init_time(obj)
            init_time = obj.init_time;
        end

        %% Set the initial time string.
        % Inputs:
        % init_time : String, formatted string representing the scenario
        %             start time, follows format_in.
        function set.init_time(obj, init_time)
            obj.init_time = init_time;
            obj.init_time_datetime = datetime(obj.init_time, ...
                'InputFormat', obj.format_in, 'TimeZone', 'UTC');
        end

        %% Convert input string into datetime object
        % Inputs:
        % time      : String, formatted string representing a time value,
        %             formatted to follow format_in.
        % Outputs:
        % dt_obj    : Datetime object, converted time input
        function dt_obj = time_2_datetime(obj, time)
            dt_obj = datetime(time, 'InputFormat', obj.format_in,...
                'TimeZone', 'UTC');
        end

        %% Convert formatted time string to EpSec.
        % Inputs:
        % time  : String, formatted time string following format_in
        % Outputs:
        % diff  : Float, difference between scenario start time (init_time)
        %         and input time string in EpSec.
        function diff = time_2_epsec(obj, time)
            diff = seconds(datetime(time, 'InputFormat', obj.format_in, ...
                'TimeZone', 'UTC') - obj.init_time_datetime);
        end

        %% Convert value in EpSec into datetime object.
        % Inputs:
        % time      : Float, time in EpSec (seconds since scenario start
        %             time)
        % Outputs:
        % time_new  : datetime object, time as a datetime object.
        function time_new = epsec_2_time(obj, time)
            time_new = obj.init_time_datetime + duration(0, 0, time);
        end

        %% Converts value in EpSec into date vector for use with the
        % eci2lla function.
        % Inputs:
        % time      : Float, time in EpSec.
        % Outputs:
        % time_new  : Float vector, date vector for the new time
        function time_new = epsec_2_date_vec(obj, time)
            time_new = datevec(obj.init_time_datetime + duration(0, 0, ...
                time));
        end

        %% Convert value from EpSec to formatted string according to 
        % format_in.
        % Inputs:
        % time      : Float, time in EpSec
        % Outputs:
        % time_new  : String, formatted time string according to format_in
        function time_new = epsec_2_fstr(obj, time)
            time_new = datestr(obj.init_time_datetime + ...
                duration(0, 0, time), obj.format_out);
        end
    end
end