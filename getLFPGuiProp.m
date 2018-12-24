function props = getLFPGuiProp( )
%getLFPGuiProps does not require input
% returns a struct containing fields 
% - plotting function (options)



props.x = {'Frequency'};

props.y = {'Power'};

props.fct = {'Raw Power', ...
    'Stimulus Norm Power', ...
    'STA Power', ...
    'SFC Power'};


end

