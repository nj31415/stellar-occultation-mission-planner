%% Nicholas Jones - njones31@vt.edu
% Star
% The Star class contains data pertinent to defining a star within STK and
% the additional data needed for occultation modelling.

classdef Star < handle
    properties (SetAccess = public)
        name
        hd_num = 0;     % Default value
        ra
        dec
        ra_proper_motion
        dec_proper_motion
        parallax
    end
    
    methods
        %% Constructor for a star object. The HD and TD1 catalog numbers
        % are optional. Otherwise all properties must be defined. All
        % values are referenced to J2000 frame. Throws errors when
        % encountering data that violates conditions noted in Inputs below.
        %
        % Inputs:
        % name              : String, identifier for the star object. Must
        %                     follow STK object naming conventions.
        % ra                : Float, right ascension in degrees.
        %                     Must be between -180 and 360 degrees.
        % dec               : Float, declination in degrees. Must be
        %                     between -90 and 90 degrees.
        % ra_proper_motion  : Float, right ascension proper motion in
        %                     arcsec yr^-1. Must be between -100 and 100
        %                     arcsec yr^-1.
        % dec_proper_motion : Float, declination proper motion in arcsec
        %                     yr^-1. Must be between -100 and 100 arcsec
        %                     yr^-1.
        % parallax          : Float, parallax in arcsec. Must be between 0
        %                     and 3600 arcsec.
        % hd_num            : Int, HD catalog number. Optional
        % Outputs:
        % obj               : Star object
        function obj = Star(name, ra, dec, ra_proper_motion, ...
                dec_proper_motion, parallax, hd_num)
            if nargin > 0
                obj.name = name;
                obj.ra = ra;
                obj.dec = dec;
                obj.ra_proper_motion = ra_proper_motion;
                obj.dec_proper_motion = dec_proper_motion;
                obj.parallax = parallax;
                
                % If it was provided, set the HD catalog numbers
                if exist('hd_num', 'var') && ~isempty(hd_num)
                    obj.hd_num = hd_num;
                end
            end
        end

        %% Get the star name.
        % Outputs:
        % name  : String, the star name
        function name = get.name(obj)
            name = obj.name;
        end

        %% Set the star name.
        % Inputs:
        % name  : String, identifier of the star object. Must follow STK
        %         naming conventions.
        function set.name(obj, name)
            if ~contains(name, ' ')
                obj.name = name;
            else
                error('Star::set.name: Star names cannot have spaces');
            end
        end
        
        %% Get the star right ascension.
        % Outputs:
        % ra    : Float, the right ascension in degrees. Referenced to
        %         J2000.
        function ra = get.ra(obj)
            ra = obj.ra;
        end

        %% Set the star right ascension.
        % Inputs:
        % ra    : Float, the right ascension of the star in degrees. Must
        %         be between -180 and 360 degrees. Referenced to J2000.
        function set.ra(obj, ra)
            if ra >= -180 && ra <= 360
                obj.ra = ra;
            else
                error(['Star:set.ra: Right asecension falls outside ',...
                    'the acceptable bounds']);
            end
        end

        %% Get the star declination.
        % Outputs:
        % dec   : Float, the declination in degrees. Referenced to J2000.
        function dec = get.dec(obj)
            dec = obj.dec;
        end

        %% Set the star declination.
        % Inputs:
        % dec   : Float, the declination of the star in degrees. Must be
        %         between -90 and 90 degrees. Referenced to J2000.
        function set.dec(obj, dec)
            if dec >= -90 && dec <= 90
                obj.dec = dec;
            else
                error(['Star:set.dec: Declination falls outside the ', ...
                    'acceptable bounds']);
            end
        end

        %% Get the star right ascension proper motion.
        % Outputs:
        % ra_proper_motion  : Float, the right ascension proper motion of
        %                     the star in arcsec yr^-1. Referenced to
        %                     J2000.
        function ra_proper_motion = get.ra_proper_motion(obj)
            ra_proper_motion = obj.ra_proper_motion;
        end

        %% Set the star right ascension proper motion.
        % Inputs:
        % ra_proper_motion  : Float, the right ascension proper motion of
        %                     the star in arcsec yr^-1. Referenced to
        %                     J2000. Must be between -100 and 100 arcsec
        %                     yr^-1.
        function set.ra_proper_motion(obj, ra_proper_motion)
            if ra_proper_motion >= -100 && ra_proper_motion <= 100
                obj.ra_proper_motion = ra_proper_motion;
            else
                error(['Star::set.ra_proper_motion: Right ascension ', ...
                    'proper motion falls outside the acceptable bounds']);
            end
        end

        %% Get the star declination proper motion.
        % Outputs:
        % dec_proper_motion : Float, the declination proper motion of the
        %                     star in arcsec yr^-1. Referenced to J2000.
        function dec_proper_motion = get.dec_proper_motion(obj)
            dec_proper_motion = obj.dec_proper_motion;
        end

        %% Set the star declination proper motion.
        % Inputs:
        % dec_proper_motion : Float, the declination proper motion of the
        %                     star in arcsec yr^-1. Referenced to J2000.
        %                     Must be between -100 and 100 arc yr^-1.
        function set.dec_proper_motion(obj, dec_proper_motion)
            if dec_proper_motion >= -100 && dec_proper_motion <= 100
                obj.dec_proper_motion = dec_proper_motion;
            else
                error(['Star::set.de_proper_motion: Declination proper',...
                    ' motion falls outside the acceptable bounds']);
            end
        end

        %% Get the stellar parallax.
        % Outputs:
        % parallax  : Float, the stellar parallax in arcsec. Referenced to
        %             J2000.
        function parallax = get.parallax(obj)
            parallax = obj.parallax;
        end

        %% Set the stellar parallax.
        % InputsL
        % parallax  : Float, the stellar parallax in arcsec. Referenced to
        %             J2000. Must be between 0 and 36000 arcsec.
        function set.parallax(obj, parallax)
            if parallax >= 0 && parallax <= 3600
                obj.parallax = parallax;
            else
                error(['Star::set.parallax: Parallax falls outside the',...
                    ' acceptable bounds']);
            end
        end

        %% get the HD number.
        % Outputs:
        % hd_num    : Int, HD number of the star. A value of 0 is the
        %             default, and indicates it has not been set.
        function hd_num = get.hd_num(obj)
            hd_num = obj.hd_num;
        end

        %% Set the HD number.
        % Inputs:
        % hd_num    : Int, HD number of the star.
        function set.hd_num(obj, hd_num)
            obj.hd_num = hd_num;
        end

        %% Checks equality of two star objects. Star objects are equal if
        % name, right ascension, and declination are equal.
        %
        % Inputs:
        % obj       : The calling star object.
        % other     : Star object to compare obj to.
        %
        % Outputs:
        % tf_res    : boolean, true false value for equality of the inputs.
        function tf_res = equals(obj, other)
            tf_res = false;
            if exist('obj', 'var') && exist('other', 'var') && ...
                strcmp(class(obj), class(other))
                tf_res = obj == other || strcmp(obj.name, other.name) &&...
                    obj.ra == other.ra && obj.dec == other.dec;
            end
        end
    end
end