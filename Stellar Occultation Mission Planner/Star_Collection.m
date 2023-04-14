%% Nicholas Jones - njones31@vt.edu
% Star_Collection
% The star collection class contains data pertinent to defining stars
% within STK and additional data needed for occultation modelling. It also
% enables seaching stars by location in right ascension and declination.
% Utilizes the J2000 reference frame.

classdef Star_Collection < handle
    properties (SetAccess = private)
        % List of star objects
        star_list = Star.empty;
    end
    
    methods
        %% Constructor for Star_Collection object. All stellar coordinates
        % are referenced to J2000.
        % name_vector   : String array, contains star name(s)
        % data_matrix   : Float matrix, contains data for stars. Should be
        %                 formatted by column as:
        %                   1. HD catalog number
        %                   2. Right ascension - degrees
        %                   3. Declination - degrees
        %                   4. Right ascension proper motion - arcsec yr^-1
        %                   5. Declination proper motion - arcsec yr^-1
        %                   6. Parallax - arcsec
        function obj = Star_Collection(name_vector, data_matrix)
            obj.add_star(name_vector, data_matrix);
        end

        %% Get the names of all stars in the collection.
        % Outputs:
        % name_vector   : String array, contains all the star names
        function name_vector = get_name_vector(obj)
            name_vector = [obj.star_list.name];
        end

        %% Get the star_data_matrix for all stars in the collection.
        % Outputs:
        % data_matrix   : Float vector, contains data for all stars in the
        % star collection
        function data_matrix = get_data(obj)
            data_matrix = [obj.star_list.hd_num obj.star_list.ra ...
                obj.star_list.dec obj.star_list.ra_proper_motion...
                obj.star_list.dec_proper_motion obj.star_list.parallax];
        end
        
        %% Add a star to the collection
        % Inputs:
        % name  : String array, name(s) of the star(s) to add to the 
        %         collection. If inputting a string array, should be a
        %         column.
        % data  : Float matrix, star data formatted consistently with
        %         constructor format
        function add_star(obj, name, data)
            if size(name, 1) ~= size(data, 1)
                error(['Star_Collection::add_star: Length of names does'...
                    'not match data']);
            elseif size(data, 2) ~= 6
                error(['Star_Collection::add_star: Incorrect data'...
                    'formatting']);
            end
            
            for i = 1 : size(name, 1)
                if isempty(obj.star_list) || ...
                        ~any(matches([obj.star_list.name], name(i)))
                    obj.star_list(end + 1, 1) = Star(name(i), ...
                        data(i, 2), data(i, 3), data(i, 4), data(i, 5), ...
                        data(i, 6), data(i, 1));
                else
                    error(['Star_Collection::add_star: Duplicate star', ...
                        ' name detected. Faield at : ' num2str(i)]);
                end
            end
        end

        %% Remove a star from the collection.
        % Inputs:
        % name  : String array, name(s) of the star(s) to remove from the
        %         collection. If inputting a string array, should be a
        %         column.
        function remove_star(obj, name)
            remove_idx = matches([obj.star_list.name], name);
            obj.star_list(remove_idx) = [];
        end

        %% Search the star collection by star names.
        % Inputs:
        % name          : String array, name(s) of the star(s) to search
        %                 the collection for.
        % Outputs:
        % star_list     : Star array, contains found stars
        function star_list = search_names(obj, name)
            search_idx = matches([obj.star_list.name], name);
            star_list = obj.star_list(search_idx);
        end

        %% Search the star collection by right ascension and declination.
        % Search limits are inclusive.
        % Inputs:
        % ra_lim        : Float vector, specification of the min and max
        %                 right ascension in deg. [min_ra max_ra]
        % dec_lim       : Float vector, specification of the min and max
        %                 declination in deg. [min_dec max_dec]
        % Outputs
        % star_list     : Star array, contains found stars
        function star_list = ...
                search_region(obj, ra_lim, dec_lim)
            search_idx = [obj.star_list.ra] >= ra_lim(1) &...
                [obj.star_list.ra] <= ra_lim(2) &...
                [obj.star_list.dec] >= dec_lim(1) &...
                [obj.star_list.dec] <= dec_lim(2);

            star_list = obj.star_list(search_idx);
        end

        %% Function to get the number of stars in the collection
        % Outputs:
        % num_star  : Int, number of stars in the collection
        function num_star = count(obj)
            num_star = length(obj.star_list);
        end

        %% Function to get a star at an index in the collection
        % Inputs:
        % idx   : Int, star index
        % Outputs:
        % star  : Star, star at requested index
        function star = get_star(obj, idx)
            star = obj.star_list(idx);
        end
    end
end