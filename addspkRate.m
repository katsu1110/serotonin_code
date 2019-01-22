
%
% add...
% - ...number of spikes
% - ...number of presented frames
% - ...spike rate
%
% @author Corinna Lorenz
% created 22.09.2015
% changed 22.09.2015
% adapted 11.04.2017 - added the option to compute the spikes/cycle/seconds
%


for n = 1:length(ex.Trials)
    
    % only consider correctly executed trials
    if ex.Trials(n).Reward == 1
        t_spks = ex.Trials(n).Spikes;   % time of spike occurance
        f_strt = ex.Trials(n).Start;    % relative time of frame onset
        t_strt = f_strt - ex.Trials(n).TrialStart; % aligned time of stimulus frames onset
        
        % if it is an experiment with adaptation, you should only consider the
        % time after the neuron adapted to the stimulus, i.e. after the
        % initiated adaption duration
        if ex.stim.vals.adaptation
            ex.Trials(n).adapt = true;
            t_strt = t_strt(t_strt > ex.stim.vals.adaptationDur);
        else
            ex.Trials(n).adapt = false;
            ex.Trials(n).adaptationDuration = 0;
        end
        
        
        frame_dur = mean(diff(t_strt)); % average frame duration
        t_end = t_strt(end)+frame_dur;  % time of stimulus ending
        
        ex.Trials(n).stimDuration   = t_end - t_strt(1); % stimulus duration
        ex.Trials(n).spkCount       = length(find( t_spks>=t_strt(1) & t_spks<=t_end)); % spikes occuring during stimulus presentation
        ex.Trials(n).spkRate        = ex.Trials(n).spkCount / (t_strt(end) - t_strt(1));  % resulting spike rate
        
        
        %% only uncomment if you want the spike rate computed over a full stimulus cycle
        % note: this reduces the analysis window.
        
        %     % Here, we consider the response per cycle/s. Since simple cells are
        %     % sensitive to the phase of the grating and it is more complicated to
        %     % compute the temporal response pattern, we decided to only
        %     % consider the stimulus presentation of full cycles.
        %     % Also keep in mind that the starting phase was not randomized for
        %     % experiments with drifting gratings. Therefore, we have to keep in
        %     % mind the response latency that confounds the spike counts in the
        %     % beginning of the trial. To overcome this problem we assume a lag of
        %     % 50ms and start the analysis there.
        %
        %     resplatency = 0.05; % we estimate the response latency to be 50ms
        %     tf = ex.stim.vals.tf; % temporal frequency of the drifting grating
        %
        %     n_cycles = floor( ex.Trials(n).stimDuration * tf ); % number of full cycles
        %
        %     if tf>2.5
        %         t_strt = t_strt+resplatency;
        %         t_end = t_strt(1) + n_cycles/tf; % start and end of the response matching n full cycles
        %         ex.Trials(n).spkCountAllFullCycles = length(find( t_spks>=t_strt(1) & t_spks<=t_end)) ;
        %         ex.Trials(n).spkRate = ex.Trials(n).spkCountAllFullCycles / n_cycles / (1/tf)  ;
        %     end
        
        
    else
        
        ex.Trials(n).stimDuration   = -1;
        ex.Trials(n).spkCount       = -1;
        ex.Trials(n).spkRate        = -1;
        
    end
    
    
end
% for debugging
% fprintf('tf %1.0f ==> number of full cycles = %1.0f \n', tf, n_cycles);

clearvars t_spks t_strt t_end n frame_dur