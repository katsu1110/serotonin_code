function h = LFPGui( dat )
%lfpgui is the gui function to plot different lfp depending measures
% the function is called from mugui that directly opens a new window to use
% the selected data and analyse their lfp changes
%
% The aim is to be able to plot:
% - average power spectrum 5HT vs baseline
% - spike triggered average of LFP
% - spike field coherence 
% all for raw and stimulus normalized data



% preamble
guiprop = getLFPGuiProp();


%% figure
sz = get(0, 'Screensize');
fig_h = figure('Position', [sz(3)*0.2 sz(4)*0.3 sz(3)*0.4 sz(4)*0.4], ...
    'HandleVisibility', 'callback');
pos = get(fig_h, 'position');



%------------------------- popmenu for function
popFct_h = uicontrol(fig_h, ....
    'Style',  'popupmenu',...
    'String',  guiprop.fct,...
    'Position', [pos(3)*0.8 pos(4)*0.8 pos(3)*0.15 pos(4)*0.1],...
    'Callback', '');

%------------------------- popmenu for x
popX_h = uicontrol(fig_h, ....
    'Style',  'popupmenu',...
    'String', [{'x'} guiprop.x],...
    'Position', [pos(3)*0.8 pos(4)*0.7 pos(3)*0.15 pos(4)*0.1],...
    'Callback', '');

%-------------------------- popmenu for y
popY_h = uicontrol(fig_h, ....
    'Style',  'popupmenu',...
    'String', [{'y'} guiprop.y],...
    'Position', [pos(3)*0.8 pos(4)*0.6 pos(3)*0.15 pos(4)*0.1],...
    'Callback', '');

%---------------------------- plot button
but_plot = uicontrol(fig_h,...
    'Style', 'pushbutton',...
    'Position', [pos(3)*0.8 pos(4)*0.1 pos(3)*0.15 pos(4)*0.05],...
    'String',   'plot', ...
    'Callback', @callPlot);

%----------------------------checkbox for different stimulus types
popStimulus_h = uicontrol(fig_h, ....
    'Style',  'checkbox',...
    'String', 'Stimulus Dep',...
    'Position', [pos(3)*0.8 pos(4)*0.4 pos(3)*0.15 pos(4)*0.1],...
    'Callback', '');


%% 

    function callPlot(~, ~, ~)
        
        % whether for different stimulus parameters (than color coded
        % power) or normal spectrum         
        stimnorm = get(popStimulus_h, 'Value');
        
                
        switch guiprop.fct{get(popFct_h, 'Value')}
            case 'Raw Power'
                rawpow(dat, 'stimnorm', stimnorm);
            case 'Stimulus Norm Power'
                stimpow(dat, 'stimnorm', stimnorm);
            case 'STA Power'
                spktrig(dat, 'stimnorm', stimnorm);
            case 'SFC Power'
                spkfieldco(dat, 'stimnorm', stimnorm);
        end
        
        
        set(gca, 'Position', [.15 .1 .5 .8]);
        axis square
        box off
        
    end



end






