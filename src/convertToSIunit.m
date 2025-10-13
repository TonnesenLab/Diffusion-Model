function scale=convertToSIunit(unit)
    switch unit
        case 's'
            scale = 1;
        case 'ms'
            scale = 10^-3;
        case 'µs'
            scale = 10^-6;
        case 'M'
            scale = 10^3;
        case 'mM'
            scale = 1;
        case 'µM'
            scale = 10^-3;
        case 'nM'
            scale = 10^-6;
    end
end