function ex = filterLFP(ex)
% filter LFP traces into 'band'

% sampling rate
fs = 1000;
fsignal = 'MPsignal';

% filters
bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
range = {[0.5,4], [4, 7], [8, 13], [14, 29], [30, 80]};
lenb = length(bands);
filtdim = 2;
for i = 1:lenb
    [filters.(bands{i}).b, filters.(bands{i}).a] = butter(filtdim, range{i}/(fs/2), 'bandpass');
end

% assign each trial
for i = 1:length(ex.Trials)    
    for b = 1:lenb
        % all
        ex.Trials(i).(['lfp_' bands{b} '_tc']) = ...
                filter(filters.(bands{b}).b, filters.(bands{b}).a, ex.Trials(i).(fsignal));
        ex.Trials(i).(['lfp_' bands{b} '_pow']) = ...
            bandpower(ex.Trials(i).(fsignal), fs, range{b});
        % baseline, stimulus evoked, sustained
        for u = 1:3
            ex.Trials(i).period(u).(['lfp_' bands{b} '_tc']) = ...
                filter(filters.(bands{b}).b, filters.(bands{b}).a, ex.Trials(i).period(u).(fsignal));   
            ex.Trials(i).period(u).(['lfp_' bands{b} '_pow']) = ...
                bandpower(ex.Trials(i).period(u).(fsignal), fs, range{b});
        end
    end
end
