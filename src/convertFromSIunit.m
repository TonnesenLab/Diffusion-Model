function scale=convertFromSIunit(unit)
% convertFromSIunit - This function convert from the international system
% units (Time and concentraion) into the desired units. 
%Input: unit - the desired unit specified as: 's', 'ms', 'µs' for time or
%               'M', 'mM','µM','nM' for concentration
%Output: scale - a value by wich the magnitud needs to be multiplied to
%               convert it to the desired units

    switch unit
        case 's'
            scale = 1;
        case 'ms'
            scale = 10^3;
        case 'µs'
            scale = 10^6;
        case 'M'
            scale = 10^-3;
        case 'mM'
            scale = 1;
        case 'µM'
            scale = 10^3;
        case 'nM'
            scale = 10^6;
    end
end