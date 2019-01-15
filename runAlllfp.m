function runAlllfp(type)

if nargin < 1; type = 'all'; end

% path =======================================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/serotonin_project']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

% LFP filtering =============================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'filter'))   
    % load data and lists
    a = load([mypath '/Katsuhisa/serotonin_project/dataset/Data/exinfo.mat'], 'exinfo');
    exinfo = a.exinfo;
    lene = length(exinfo);
    lfplist = cell(1, lene);
    savepath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    parfor i = 1:lene
        % filenames
        if strcmp('Z:\data\mango\0132\ma_0132_c1_sortHN_all.grating.ORxRC.mat', exinfo(i).fname)
            exinfo(i).fname = 'Z:\data\mango\0132\ma_0132_c1_sortHN_4.25PM.grating.ORxRC.mat';
            exinfo(i).fname_drug = 'Z:\data\mango\0132\ma_0132_c1_sortHN_4.30PM.grating.ORxRC_5HT.mat';
        end
        
        slash = strfind(exinfo(i).fname, '\');
        fname = exinfo(i).fname(slash(end)+1:end);        
        slash = strfind(exinfo(i).fname_drug, '\');
        fname_drug = exinfo(i).fname_drug(slash(end)+1:end);        
        
        lfplist{i} = {fname, fname_drug};
        
        % baseline 
        try
            loadCluster(exinfo(i).fname, 'loadlfp', 1, ...
                'save', [savepath fname]);
        catch
            lfplist{i} = {nan, nan};
            disp([fname ': error in loadCluster for baseline. skip...'])
            continue
        end        
        
        % drug
        try
            loadCluster(exinfo(i).fname_drug, 'loadlfp', 1, ...
                'save', [savepath fname_drug]);
        catch
            lfplist{i} = {nan, nan};
            disp([fname_drug ': error in loadCluster for drug. skip...'])
            continue
        end
    end
    save([savepath 'lfplist.mat'], 'lfplist')
    disp('all saved!')
end


% MP preprocessing (stimulus response) =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'MP'))
    % path
    addpath(genpath([mypath '/Katsuhisa/code/integrated/MP']))
    
    % load data and lists
    a = load([mypath '/Katsuhisa/serotonin_project/dataset/Data/exinfo.mat'], 'exinfo');
    exinfo = a.exinfo;
    lene = length(exinfo);
    N = lene;
    for i = 1:N   
        % filenames
        if strcmp('Z:\data\mango\0132\ma_0132_c1_sortHN_all.grating.ORxRC.mat', exinfo(i).fname)
            exinfo(i).fname = 'Z:\data\mango\0132\ma_0132_c1_sortHN_4.25PM.grating.ORxRC.mat';
            exinfo(i).fname_drug = 'Z:\data\mango\0132\ma_0132_c1_sortHN_4.30PM.grating.ORxRC_5HT.mat';
        end
        slash = strfind(exinfo(i).fname, '\');
        fname = exinfo(i).fname(slash(end)+1:end);   
        slash = strfind(exinfo(i).fname_drug, '\');
        fname_drug = exinfo(i).fname_drug(slash(end)+1:end);        
        
        if contains(fname_drug, 'xRC') % RC only
%             try
                % baseline 
                ex0 = loadCluster(exinfo(i).fname, 'loadlfp', 1, ...
                    'filtering', 0);        

                % drug
                ex2 = loadCluster(exinfo(i).fname_drug, 'loadlfp', 1, ...
                     'filtering', 0);
            
                % perform MP & save the ex
                MP_single(ex0, 0, fname);
                MP_single(ex2, 0, fname_drug);
%             catch
%                 disp(['Err: ' fname_drug])
%             end
        end
    end
end

% MP preprocessing (internally generated LFP) =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'MP1'))
    % path
    addpath(genpath([mypath '/Katsuhisa/code/integrated/MP']))
    
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');
    lfplist = a.lfplist;
    N = length(lfplist);
    for i = 1:N
        try
            % load data
            base = load([loadpath lfplist{i}{1}], 'ex');
            drug = load([loadpath lfplist{i}{2}], 'ex');

            % perform MP & save the ex
            MP_single(base.ex, 1, lfplist{i}{1});
            MP_single(drug.ex, 1, lfplist{i}{2});
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
end


% MP structure =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'MPs'))    
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
    
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    loadpath_mp = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/MPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
     % initialization
    Out1 = cell(1, N); Out2 = cell(1, N);
    Out3 = cell(1, N); Out4 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            goodunit(i) = ismember(i, incl_i);
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            % load =======================            
            % load lfp data
            base = load([loadpath lfplist{i}{1}], 'ex');
            drug = load([loadpath lfplist{i}{2}], 'ex');
            
            % load MP data
            mp_base = load([loadpath_mp 'Trials/' lfplist{i}{1}], 'exn');
            mp_drug = load([loadpath_mp 'Trials/' lfplist{i}{2}], 'exn');
            mp_base1 = load([loadpath_mp 'iTrials/' lfplist{i}{1}], 'exn');
            mp_drug1 = load([loadpath_mp 'iTrials/' lfplist{i}{2}], 'exn');
            
            % analysis
            Out1{i} = MPspectrogram(base.ex, mp_base.exn);
            Out2{i} = MPspectrogram(drug.ex, mp_drug.exn);
            Out3{i} = MPspectrogram(base.ex, mp_base1.exn);
            Out4{i} = MPspectrogram(drug.ex, mp_drug1.exn);
            
            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, base.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    
    % unit info
    oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
    Mps_all.lfplist = lfplist(oks);
    Mps_all.stmtype = stmtype(oks);
    Mps_all.animal = animal(oks);
    Mps_all.is5ht = is5ht(oks);
    Mps_all.goodunit = goodunit(oks);

    % results from analysis    
    iMps_all = Mps_all;
    Mps_all.base = Out1(oks); 
    Mps_all.drug = Out2(oks); 
    iMps_all.base = Out3(oks); 
    iMps_all.drug = Out4(oks); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit', 'base', 'drug'};
    for s = 1:length(stmtypes)
        idx = Mps_all.stmtype == s;
        
        for f = 1:length(fields)
            Mps.(fields{f}) = Mps_all.(fields{f})(idx);
            iMps.(fields{f}) = iMps_all.(fields{f})(idx);
        end
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/MPs/Mps_' stmtypes{s} '.mat'], 'Mps', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/MPs/iMps_' stmtypes{s} '.mat'], 'iMps', '-v7.3')
        disp([stmtypes{s} ': Mps saved!'])
    end
end


% IRASA preprocessing =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'IRASA'))
    % path
    addpath(genpath([mypath '/Katsuhisa/code/integrated/IRASA']))
    
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');
    lfplist = a.lfplist;
    N = length(lfplist);
    parfor i = 1:N
        try
            if contains(lfplist{i}{1}, 'xRC')
                % load data
                base = load([loadpath lfplist{i}{1}], 'ex');
                drug = load([loadpath lfplist{i}{2}], 'ex');
                
                % z-scoring
                [ex0, ex2] = zscore_pairLFP(base.ex, drug.ex, 'LFP_prepro');
                [ex0, ex2] = zscore_pairLFP(ex0, ex2, 'iLFP_prepro');
                
                % perform IRASA & save the ex
                IRASA_single(ex0, lfplist{i}{1});
                IRASA_single(ex2, lfplist{i}{2});
            end
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
end

% % basic LFP analysis =========================
% if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps'))
%     % stimulus types
%     stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
%    
%     % load data and lists
%     loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
%     a = load([loadpath 'lfplist.mat'], 'lfplist');  
%     b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
%     c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
%     incl_i = c.incl_i;
%     list_RC = b.list_RC;
%     lfplist = a.lfplist;
%     N = length(lfplist);
%     
%     % initialization
%     Out1 = cell(1, N); Out2 = cell(1, N);
%     Out3 = cell(1, N); Out4 = cell(1, N);
%     goodunit = zeros(1, N); is5ht = zeros(1, N);
%     stmtype = zeros(1, N); animal = zeros(1, N); 
%     parfor i = 1:N            
%         % goodunit
%         if isempty(strfind(lfplist{i}{1}, 'xRC'))
%             goodunit(i) = ismember(i, incl_i);
%         else
%             goodunit(i) = ismember(i, list_RC);
%             stmtype(i) = 1;
%         end
%         try
%             % baseline =======================
%             d0 = load([loadpath lfplist{i}{1}], 'ex');
%             Out1{i} = stmLFP(d0.ex, 'LFP_prepro', {'all'});
%             Out2{i} = stmLFP(d0.ex, 'iLFP_prepro', {'all'});
% 
%             % stimlus type
%             if stmtype(i) == 0
%                 stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
%             end
% 
%             % is mango
%             animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');
% 
%             % drug ===========================
%             d2 = load([loadpath lfplist{i}{2}], 'ex');
%             Out3{i} = stmLFP(d2.ex, 'LFP_prepro', {'all'});
%             Out4{i} = stmLFP(d2.ex, 'iLFP_prepro', {'all'});
% 
%             % is 5HT
%             is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));
% 
%             disp(['session ' num2str(i) ' analyzed!'])
%         catch
%             disp(['session ' num2str(i) ' error'])
%         end
%     end
%     % unit info
%     lfplist = lfplist(~cellfun('isempty', Out1));
%     Lfps_all.lfplist = lfplist;
%     Lfps_all.stmtype = stmtype(~cellfun('isempty', Out1));
%     Lfps_all.animal = animal(~cellfun('isempty', Out1));
%     Lfps_all.is5ht = is5ht(~cellfun('isempty', Out1));
%     Lfps_all.goodunit = goodunit(~cellfun('isempty', Out1));
% 
%     % results from analysis    
%     iLfps_all = Lfps_all;
%     Lfps_all.cond(1).LFP_prepro = Out1(~cellfun('isempty', Out1)); 
%     Lfps_all.cond(2).LFP_prepro = Out3(~cellfun('isempty', Out3)); 
%     iLfps_all.cond(1).LFP_prepro = Out2(~cellfun('isempty', Out2)); 
%     iLfps_all.cond(2).LFP_prepro = Out4(~cellfun('isempty', Out4)); 
%     
%     % split by stimulus type
%     fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
%     for s = 1:length(stmtypes)
%         idx = Lfps_all.stmtype == s;
%         
%         for f = 1:length(fields)
%             Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
%             iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
%         end
%         
%         for d = 1:2
%             Lfps.cond(d).LFP_prepro = Lfps_all.cond(d).LFP_prepro(idx);
%             iLfps.cond(d).LFP_prepro = iLfps_all.cond(d).LFP_prepro(idx);
%         end
%         
%         % autosave
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps/Lfps_' stmtypes{s} '.mat'], 'Lfps', '-v7.3')
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps/iLfps_' stmtypes{s} '.mat'], 'iLfps', '-v7.3')
%         disp([stmtypes{s} ': lfps saved!'])
%     end
% end

% synchronization index =====================
if sum(strcmp(type, 'all')) || sum(strcmp(type, 'SI'))
    % load data and lists
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    a = load([mypath '/Katsuhisa/serotonin_project/dataset/Data/exinfo.mat'], 'exinfo');
    exinfo = a.exinfo;
    N = length(exinfo);        
        
    % initialization
    Out1 = cell(1, N); Out2 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % filenames
        fname = strrep(exinfo(i).fname, '\', '/');
        fname = [mypath fname(3:end)];
        fname_drug = strrep(exinfo(i).fname_drug, '\', '/');
        fname_drug = [mypath fname_drug(3:end)];
        
        % goodunit
        if contains(fname, 'xRC')
            goodunit(i) = ismember(i, incl_i);
            analysiswnd = [0.8 0];
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
            analysiswnd = [0.2 0];
        end
        try
            % analysis
            Out1{i} = getSI(fname, analysiswnd, 50, 1000);
            Out2{i} = getSI(fname_drug, analysiswnd, 50, 1000);
            
            % stimlus type
            if stmtype(i) == 0
                if ismember(1, contains(fname, '.OR'))
                    stmtype(i) = 2;
                elseif ismember(1, contains(fname, '.CO'))
                    stmtype(i) = 3;
                elseif ismember(1, contains(fname, '.SF'))
                    stmtype(i) = 4;
                elseif ismember(1, contains(fname, '.SZ'))
                    stmtype(i) = 5;  
                end
            end

            % is mango
            animal(i) = contains(fname, 'mango');

            % is 5HT
            is5ht(i) = ismember(1, contains(fname_drug, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
%     SI_all.fname{1} = fname(oks);
%     SI_all.fname{2} = fname_drug(oks);
    SI_all.stmtype = stmtype(oks);
    SI_all.animal = animal(oks);
    SI_all.is5ht = is5ht(oks);
    SI_all.goodunit = goodunit(oks);

    % results from analysis    
    SI_all.si{1} = Out1(oks); 
    SI_all.si{2} = Out2(oks);
            
    % autosave
    save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/SI/SI_all.mat'], 'SI_all', '-v7.3')
    disp('SI_all saved!')
end

% basic spike count analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'spike_count'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            goodunit(i) = ismember(i, incl_i);
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            % load =======================
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            d2 = load([loadpath lfplist{i}{2}], 'ex');
            
            % analysis
            Out1{i} = pair_stmLFP(d0.ex, d2.ex, 'LFP_prepro', 0, {'mat'}, 1);

            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    oks = ~cellfun('isempty', Out1);
    trmat.lfplist = lfplist(oks);
    trmat.stmtype = stmtype(oks);
    trmat.animal = animal(oks);
    trmat.is5ht = is5ht(oks);
    trmat.goodunit = goodunit(oks);

    % results from analysis    
    trmat.LFP_prepro = Out1(oks);  
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = trmat.stmtype == s;
        
%         if s > 1
%             continue
%         end
        
        for f = 1:length(fields)
            Trmat.(fields{f}) = trmat.(fields{f})(idx);
        end
        
        Trmat.LFP_prepro = trmat.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/trialmat/Trmat_' stmtypes{s} '.mat'], 'Trmat', '-v7.3')
        disp([stmtypes{s} ': Trmat saved!'])
    end
end

% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); Out2 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            continue
%             goodunit(i) = ismember(i, incl_i);
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            % load =======================
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            d2 = load([loadpath lfplist{i}{2}], 'ex');
            
            % analysis
            Out1{i} = pair_stmLFP(d0.ex, d2.ex, 'LFP_prepro');
            Out2{i} = pair_stmLFP(d0.ex, d2.ex, 'iLFP_prepro');

            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
    Lfps_all.lfplist = lfplist(oks);
    Lfps_all.stmtype = stmtype(oks);
    Lfps_all.animal = animal(oks);
    Lfps_all.is5ht = is5ht(oks);
    Lfps_all.goodunit = goodunit(oks);

    % results from analysis    
    iLfps_all = Lfps_all;
    Lfps_all.LFP_prepro = Out1(oks); 
    iLfps_all.LFP_prepro = Out2(oks); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = Lfps_all.stmtype == s;
        
        if s > 1
            continue
        end
        
        for f = 1:length(fields)
            Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
            iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
        end
        
        Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
        iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair/Lfps_' stmtypes{s} '.mat'], 'Lfps', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair/iLfps_' stmtypes{s} '.mat'], 'iLfps', '-v7.3')
        disp([stmtypes{s} ': lfps saved!'])
    end
end

% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_nothin'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); Out2 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            goodunit(i) = ismember(i, incl_i);
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            % load =======================
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            d2 = load([loadpath lfplist{i}{2}], 'ex');
            
            % analysis (no thinning)
            Out1{i} = pair_stmLFP(d0.ex, d2.ex, 'LFP_prepro', 0);
            Out2{i} = pair_stmLFP(d0.ex, d2.ex, 'iLFP_prepro', 0);

            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info    
    oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
    Lfps_all.lfplist = lfplist(oks);
    Lfps_all.stmtype = stmtype(oks);
    Lfps_all.animal = animal(oks);
    Lfps_all.is5ht = is5ht(oks);
    Lfps_all.goodunit = goodunit(oks);

    % results from analysis    
    iLfps_all = Lfps_all;
    Lfps_all.LFP_prepro = Out1(oks); 
    iLfps_all.LFP_prepro = Out2(oks); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = Lfps_all.stmtype == s;
        
        for f = 1:length(fields)
            Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
            iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
        end
        
        Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
        iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin/Lfps_' stmtypes{s} '.mat'], 'Lfps', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin/iLfps_' stmtypes{s} '.mat'], 'iLfps', '-v7.3')
        disp([stmtypes{s} ': lfps saved!'])
    end
end


% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_ps'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); Out2 = cell(1, N);
    Out3 = cell(1, N); Out4 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            continue
%             goodunit(i) = ismember(i, incl_i);
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            % load =======================
            % baseline ---------------------------
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            
            % split trials based on pupil
            ex = ex_spliter(d0.ex);
            
            % analysis
            Out1{i} = pair_stmLFP(ex{1}, ex{2}, 'LFP_prepro');
            Out2{i} = pair_stmLFP(ex{1}, ex{2}, 'iLFP_prepro');
            
            % drug -----------------------------------
            d2 = load([loadpath lfplist{i}{2}], 'ex');
            
            % split trials based on pupil
            ex = ex_spliter(d2.ex);
            
            % analysis
            Out3{i} = pair_stmLFP(ex{1}, ex{2}, 'LFP_prepro');
            Out4{i} = pair_stmLFP(ex{1}, ex{2}, 'iLFP_prepro');
            
            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    outs = cellfun('isempty', Out1) | cellfun('isempty', Out2) | ...
        cellfun('isempty', Out3) | cellfun('isempty', Out4);
    lfplist = lfplist(outs==0);
    Lfps_all.lfplist = lfplist;
    Lfps_all.stmtype = stmtype(outs==0);
    Lfps_all.animal = animal(outs==0);
    Lfps_all.is5ht = is5ht(outs==0);
    Lfps_all.goodunit = goodunit(outs==0);

    % results from analysis    
    iLfps_all = Lfps_all;
    Lfps_dr = Lfps_all;
    iLfps_dr = Lfps_all;
    Lfps_all.LFP_prepro = Out1(outs==0); 
    iLfps_all.LFP_prepro = Out2(outs==0); 
    Lfps_dr.LFP_prepro = Out3(outs==0); 
    iLfps_dr.LFP_prepro = Out4(outs==0); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = Lfps_all.stmtype == s;
        
        if s > 1
            continue
        end
        
        % baseline 
        for f = 1:length(fields)
            Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
            iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
        end
        
        Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
        iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_ps/Lfps_' stmtypes{s} '_base.mat'], 'Lfps', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_ps/iLfps_' stmtypes{s} '_base.mat'], 'iLfps', '-v7.3')
        
        % drug
        for f = 1:length(fields)
            Lfps.(fields{f}) = Lfps_dr.(fields{f})(idx);
            iLfps.(fields{f}) = iLfps_dr.(fields{f})(idx);
        end
        
        Lfps.LFP_prepro = Lfps_dr.LFP_prepro(idx);
        iLfps.LFP_prepro = iLfps_dr.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_ps/Lfps_' stmtypes{s} '_drug.mat'], 'Lfps', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_ps/iLfps_' stmtypes{s} '_drug.mat'], 'iLfps', '-v7.3')
        
        disp([stmtypes{s} ': lfps saved!'])
    end
end

% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_nothin_ps'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); Out2 = cell(1, N);
    Out3 = cell(1, N); Out4 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            continue
%             goodunit(i) = ismember(i, incl_i);
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            % load =======================
            % baseline ---------------------------
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            
            % drug -----------------------------------
            d2 = load([loadpath lfplist{i}{2}], 'ex');
            
            % zscoring
            [ex0, ex2] = zscore_pairLFP(d0.ex, d2.ex, 'LFP_prepro');
            [ex0i, ex2i] = zscore_pairLFP(d0.ex, d2.ex, 'iLFP_prepro');
            
            % split trials based on pupil
            ex = ex_spliter(ex0);
            exi = ex_spliter(ex0i);
            
            % analysis
            Out1{i} = pair_stmLFP(ex{1}, ex{2}, 'LFP_prepro', 0, {'all'}, 0);
            Out2{i} = pair_stmLFP(exi{1}, exi{2}, 'iLFP_prepro', 0, {'all'}, 0);
            
            % split trials based on pupil
            ex = ex_spliter(ex2);
            exi = ex_spliter(ex2i);
            
            % analysis
            Out3{i} = pair_stmLFP(ex{1}, ex{2}, 'LFP_prepro', 0, {'all'}, 0);
            Out4{i} = pair_stmLFP(exi{1}, exi{2}, 'iLFP_prepro', 0, {'all'}, 0);
            
            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    outs = cellfun('isempty', Out1) | cellfun('isempty', Out2) | ...
        cellfun('isempty', Out3) | cellfun('isempty', Out4);
    lfplist = lfplist(outs==0);
    Lfps_all.lfplist = lfplist;
    Lfps_all.stmtype = stmtype(outs==0);
    Lfps_all.animal = animal(outs==0);
    Lfps_all.is5ht = is5ht(outs==0);
    Lfps_all.goodunit = goodunit(outs==0);

    % results from analysis    
    iLfps_all = Lfps_all;
    Lfps_dr = Lfps_all;
    iLfps_dr = Lfps_all;
    Lfps_all.LFP_prepro = Out1(outs==0); 
    iLfps_all.LFP_prepro = Out2(outs==0); 
    Lfps_dr.LFP_prepro = Out3(outs==0); 
    iLfps_dr.LFP_prepro = Out4(outs==0); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = Lfps_all.stmtype == s;
        if s > 1
            continue
        end
        % baseline 
        for f = 1:length(fields)
            Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
            iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
        end
        
        Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
        iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_ps/Lfps_' stmtypes{s} '_base.mat'], 'Lfps', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_ps/iLfps_' stmtypes{s} '_base.mat'], 'iLfps', '-v7.3')
        
        % drug
        for f = 1:length(fields)
            Lfps.(fields{f}) = Lfps_dr.(fields{f})(idx);
            iLfps.(fields{f}) = iLfps_dr.(fields{f})(idx);
        end
        
        Lfps.LFP_prepro = Lfps_dr.LFP_prepro(idx);
        iLfps.LFP_prepro = iLfps_dr.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_ps/Lfps_' stmtypes{s} '_drug.mat'], 'Lfps', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_ps/iLfps_' stmtypes{s} '_drug.mat'], 'iLfps', '-v7.3')
        
        disp([stmtypes{s} ': lfps saved!'])
    end
end

% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_sc'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    es = 0.6592;
    thre = 0.01;
    Out1 = cell(1, N); Out2 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            continue
%             goodunit(i) = ismember(i, incl_i);
%             wnd = [0.2 0];
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
            wnd = [0.8 0];
        end
        try
            % load =======================
            % baseline ---------------------------
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            
            % split trials based on spike counts
            try
                [ex0, ex2] = ex_spliter_es(d0.ex, es, thre, wnd);
            catch
                % split trials based on spike counts
                ex = ex_spliter(d0.ex, 'sc');
                ex0 = ex{2};
                ex2 = ex{1};
            end
            
            % analysis
            Out1{i} = pair_stmLFP(ex0, ex2, 'LFP_prepro');
            Out2{i} = pair_stmLFP(ex0, ex2, 'iLFP_prepro');
                    
            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    outs = cellfun('isempty', Out1) | cellfun('isempty', Out2);
    lfplist = lfplist(outs==0);
    Lfps_all.lfplist = lfplist;
    Lfps_all.stmtype = stmtype(outs==0);
    Lfps_all.animal = animal(outs==0);
    Lfps_all.is5ht = is5ht(outs==0);
    Lfps_all.goodunit = goodunit(outs==0);

    % results from analysis    
    iLfps_all = Lfps_all;
    Lfps_all.LFP_prepro = Out1(outs==0); 
    iLfps_all.LFP_prepro = Out2(outs==0); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = Lfps_all.stmtype == s;

        if s > 1
            continue
        end
        
        % baseline 
        for f = 1:length(fields)
            Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
            iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
        end
        
        Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
        iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_sc/Lfps_' stmtypes{s} '_base.mat'], 'Lfps', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_sc/iLfps_' stmtypes{s} '_base.mat'], 'iLfps', '-v7.3')
        
        disp([stmtypes{s} ': lfps saved!'])
    end
end
    
% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_nothin_sc'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    es = 0.6592;
    thre = 0.01;
    Out1 = cell(1, N); Out2 = cell(1, N);
%     Out3 = cell(1, N); Out4 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            continue
%             goodunit(i) = ismember(i, incl_i);
%             wnd = [0.2 0];
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
            wnd = [0.8 0];
        end
        try
            % load =======================
            % baseline ---------------------------
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            
            % split trials based on spike counts
            try
                [ex0, ex2] = ex_spliter_es(d0.ex, es, thre, wnd);
            catch
                % split trials based on spike counts
                ex = ex_spliter(d0.ex, 'sc');
                ex0 = ex{2};
                ex2 = ex{1};
            end
            
            % analysis
            Out1{i} = pair_stmLFP(ex0, ex2, 'LFP_prepro', 0);
            Out2{i} = pair_stmLFP(ex0, ex2, 'iLFP_prepro', 0);
            
%             % drug -----------------------------------
%             d2 = load([loadpath lfplist{i}{2}], 'ex');
            
%             % split trials based on spike counts
%             try
%                 [ex0, ex2] = ex_spliter_es(d0.ex, es, thre, wnd);
%             catch
%                 % split trials based on spike counts
%                 ex = ex_spliter(d0.ex, 'sc');
%                 ex0 = ex{2};
%                 ex2 = ex{1};
%             end
%             
%             % analysis
%             Out3{i} = pair_stmLFP(ex0, ex2, 'LFP_prepro', 0);
%             Out4{i} = pair_stmLFP(ex0, ex2, 'iLFP_prepro', 0);
            
            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    outs = cellfun('isempty', Out1) | cellfun('isempty', Out2);
%     outs = cellfun('isempty', Out1) | cellfun('isempty', Out2) | ...
%         cellfun('isempty', Out3) | cellfun('isempty', Out4);
    lfplist = lfplist(outs==0);
    Lfps_all.lfplist = lfplist;
    Lfps_all.stmtype = stmtype(outs==0);
    Lfps_all.animal = animal(outs==0);
    Lfps_all.is5ht = is5ht(outs==0);
    Lfps_all.goodunit = goodunit(outs==0);

    % results from analysis    
    iLfps_all = Lfps_all;
%     Lfps_dr = Lfps_all;
%     iLfps_dr = Lfps_all;
    Lfps_all.LFP_prepro = Out1(outs==0); 
    iLfps_all.LFP_prepro = Out2(outs==0); 
%     Lfps_dr.LFP_prepro = Out3(outs==0); 
%     iLfps_dr.LFP_prepro = Out4(outs==0); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = Lfps_all.stmtype == s;
        if s > 1
            continue
        end
        % baseline 
        for f = 1:length(fields)
            Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
            iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
        end
        
        Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
        iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_sc/Lfps_' stmtypes{s} '_base.mat'], 'Lfps', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_sc/iLfps_' stmtypes{s} '_base.mat'], 'iLfps', '-v7.3')
        
%         % drug
%         for f = 1:length(fields)
%             Lfps.(fields{f}) = Lfps_dr.(fields{f})(idx);
%             iLfps.(fields{f}) = iLfps_dr.(fields{f})(idx);
%         end
%         
%         Lfps.LFP_prepro = Lfps_dr.LFP_prepro(idx);
%         iLfps.LFP_prepro = iLfps_dr.LFP_prepro(idx);
%         
%         % autosave
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_sc/Lfps_' stmtypes{s} '_drug.mat'], 'Lfps', '-v7.3')
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_sc/iLfps_' stmtypes{s} '_drug.mat'], 'iLfps', '-v7.3')
        
        disp([stmtypes{s} ': lfps saved!'])
    end
end


% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_mp'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    loadpath_mp = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/MPprepro_RC/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); 
%     Out2 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            goodunit(i) = ismember(i, incl_i);
            continue
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            % load =======================
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            d2 = load([loadpath lfplist{i}{2}], 'ex');
            d3 = load([loadpath_mp 'Trials/' lfplist{i}{1}], 'exn');
            d4 = load([loadpath_mp 'Trials/' lfplist{i}{2}], 'exn');
%             d5 = load([loadpath_mp 'iTrials/' lfplist{i}{1}], 'exn');
%             d6 = load([loadpath_mp 'iTrials/' lfplist{i}{2}], 'exn');
            
            % analysis
            Out1{i} = pair_stmLFP4mp(d0.ex, d2.ex, d3.exn, d4.exn, 'LFP_prepro');
%             Out2{i} = pair_stmLFP4mp(d0.ex, d2.ex, d5.exn, d6.exn, 'iLFP_prepro');

            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
%     oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
    oks = ~cellfun('isempty', Out1);
    Lfps_all.lfplist = lfplist(oks);
    Lfps_all.stmtype = stmtype(oks);
    Lfps_all.animal = animal(oks);
    Lfps_all.is5ht = is5ht(oks);
    Lfps_all.goodunit = goodunit(oks);

    % results from analysis    
%     iLfps_all = Lfps_all;
    Lfps_all.LFP_prepro = Out1(oks); 
%     iLfps_all.LFP_prepro = Out2(oks); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = Lfps_all.stmtype == s;
        if s > 1
            continue
        end
        for f = 1:length(fields)
            Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
%             iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
        end
        
        Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
%         iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair/Lfps_mp_' stmtypes{s} '.mat'], 'Lfps', '-v7.3')
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair/iLfps_mp_' stmtypes{s} '.mat'], 'iLfps', '-v7.3')
        disp([stmtypes{s} ': lfps saved!'])
    end
end

% % basic LFP analysis =========================
% if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_mp_nothin'))
%     % stimulus types
%     stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
%    
%     % load data and lists
%     loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
%     loadpath_mp = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/MPprepro/'];
%     a = load([loadpath 'lfplist.mat'], 'lfplist');  
%     b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
%     c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
%     incl_i = c.incl_i;
%     list_RC = b.list_RC;
%     lfplist = a.lfplist;
%     N = length(lfplist);
%     
%     % initialization
%     Out1 = cell(1, N); 
% %     Out2 = cell(1, N);
%     goodunit = zeros(1, N); is5ht = zeros(1, N);
%     stmtype = zeros(1, N); animal = zeros(1, N); 
%     parfor i = 1:N            
%         % goodunit
%         if isempty(strfind(lfplist{i}{1}, 'xRC'))
%             goodunit(i) = ismember(i, incl_i);
%             continue
%         else
%             goodunit(i) = ismember(i, list_RC);
%             stmtype(i) = 1;
%         end
%         try
%             % load =======================
%             d0 = load([loadpath lfplist{i}{1}], 'ex');
%             d2 = load([loadpath lfplist{i}{2}], 'ex');
%             d3 = load([loadpath_mp 'Trials/' lfplist{i}{1}], 'exn');
%             d4 = load([loadpath_mp 'Trials/' lfplist{i}{2}], 'exn');
% %             d5 = load([loadpath_mp 'iTrials/' lfplist{i}{1}], 'exn');
% %             d6 = load([loadpath_mp 'iTrials/' lfplist{i}{2}], 'exn');
%             
%             % analysis
%             Out1{i} = pair_stmLFP4mp(d0.ex, d2.ex, d3.exn, d4.exn, 'LFP_prepro', 0);
% %             Out2{i} = pair_stmLFP4mp(d0.ex, d2.ex, d5.exn, d6.exn, 'iLFP_prepro', 0);
% 
%             % stimlus type
%             if stmtype(i) == 0
%                 stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
%             end
% 
%             % is mango
%             animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');
% 
%             % is 5HT
%             is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));
% 
%             disp(['session ' num2str(i) ' analyzed!'])
%         catch
%             disp(['session ' num2str(i) ' error'])
%         end
%     end
%     % unit info
% %     oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
%     oks = ~cellfun('isempty', Out1);
%     
%     Lfps_all.lfplist = lfplist(oks);
%     Lfps_all.stmtype = stmtype(oks);
%     Lfps_all.animal = animal(oks);
%     Lfps_all.is5ht = is5ht(oks);
%     Lfps_all.goodunit = goodunit(oks);
% 
%     % results from analysis    
% %     iLfps_all = Lfps_all;
%     Lfps_all.LFP_prepro = Out1(oks); 
% %     iLfps_all.LFP_prepro = Out2(oks); 
%     
%     % split by stimulus type
%     fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
%     for s = 1:length(stmtypes)
%         idx = Lfps_all.stmtype == s;
%         if s > 1
%             continue
%         end
%         for f = 1:length(fields)
%             Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
% %             iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
%         end
%         
%         Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
% %         iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
%         
%         % autosave
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin/Lfps_mp_' stmtypes{s} '.mat'], 'Lfps', '-v7.3')
% %         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin/iLfps_mp_' stmtypes{s} '.mat'], 'iLfps', '-v7.3')
%         disp([stmtypes{s} ': lfps saved!'])
%     end
% end


% % basic LFP analysis =========================
% if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_mp_ps'))
%     % stimulus types
%     stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
%    
%     % load data and lists
%     loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];    
%     loadpath_mp = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/MPprepro/'];
%     a = load([loadpath 'lfplist.mat'], 'lfplist');  
%     b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
%     c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
%     incl_i = c.incl_i;
%     list_RC = b.list_RC;
%     lfplist = a.lfplist;
%     N = length(lfplist);
%     
%     % initialization
%     Out1 = cell(1, N); 
% %     Out2 = cell(1, N);
% %     Out3 = cell(1, N);
% %     Out4 = cell(1, N);
%     goodunit = zeros(1, N); is5ht = zeros(1, N);
%     stmtype = zeros(1, N); animal = zeros(1, N); 
%     parfor i = 1:N            
%         % goodunit
%         if isempty(strfind(lfplist{i}{1}, 'xRC'))
%             goodunit(i) = ismember(i, incl_i);
%             continue
%         else
%             goodunit(i) = ismember(i, list_RC);
%             stmtype(i) = 1;
%         end
%         try
%             % load =======================
%             % baseline ---------------------------
%             d0 = load([loadpath lfplist{i}{1}], 'ex');
%             d3 = load([loadpath_mp 'Trials/' lfplist{i}{1}], 'exn');
% %             d5 = load([loadpath_mp 'iTrials/' lfplist{i}{1}], 'exn');
%             
%             % split trials based on pupil
%             [ex, idx] = ex_spliter(d0.ex);
%             exn3_0 = d3.exn;
%             exn3_0.Trials = exn3_0.Trials(idx{1});
%             exn3_2 = d3.exn;
%             exn3_2.Trials = exn3_2.Trials(idx{2});
% %             exn5_0 = d5.exn;
% %             exn5_0.Trials = exn5_0.Trials(idx{1});
% %             exn5_2 = d5.exn;
% %             exn5_2.Trials = exn5_2.Trials(idx{2});
%             
%             % analysis
%             Out1{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn3_0, exn3_2, 'LFP_prepro');
% %             Out2{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn5_0, exn5_2, 'iLFP_prepro');
%             
%             % drug -----------------------------------
% %             d2 = load([loadpath lfplist{i}{2}], 'ex');
% %             d4 = load([loadpath_mp 'Trials/' lfplist{i}{2}], 'exn');
% %             d6 = load([loadpath_mp 'iTrials/' lfplist{i}{2}], 'exn');
%             
%             % split trials based on pupil
% %             [ex, idx] = ex_spliter(d2.ex);
% %             exn4_0 = d4.exn;
% %             exn4_0.Trials = exn4_0.Trials(idx{1});
% %             exn4_2 = d4.exn;
% %             exn4_2.Trials = exn4_2.Trials(idx{2});
% %             exn6_0 = d6.exn;
% %             exn6_0.Trials = exn6_0.Trials(idx{1});
% %             exn6_2 = d6.exn;
% %             exn6_2.Trials = exn6_2.Trials(idx{2});
%             
%             % analysis
% %             Out3{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn4_0, exn4_2, 'LFP_prepro');
% %             Out4{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn6_0, exn6_2, 'iLFP_prepro');
%             
%             % stimlus type
%             if stmtype(i) == 0
%                 stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
%             end
% 
%             % is mango
%             animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');
% 
%             % is 5HT
%             is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));
% 
%             disp(['session ' num2str(i) ' analyzed!'])
%         catch
%             disp(['session ' num2str(i) ' error'])
%         end
%     end
%     % unit info
% %     outs = cellfun('isempty', Out1) | cellfun('isempty', Out2) | ...
% %         cellfun('isempty', Out3) | cellfun('isempty', Out4);
%     outs = cellfun('isempty', Out1) ;
%     lfplist = lfplist(outs==0);
%     Lfps_all.lfplist = lfplist;
%     Lfps_all.stmtype = stmtype(outs==0);
%     Lfps_all.animal = animal(outs==0);
%     Lfps_all.is5ht = is5ht(outs==0);
%     Lfps_all.goodunit = goodunit(outs==0);
% 
%     % results from analysis    
% %     iLfps_all = Lfps_all;
% %     Lfps_dr = Lfps_all;
% %     iLfps_dr = Lfps_all;
%     Lfps_all.LFP_prepro = Out1(outs==0); 
% %     iLfps_all.LFP_prepro = Out2(outs==0); 
% %     Lfps_dr.LFP_prepro = Out3(outs==0); 
% %     iLfps_dr.LFP_prepro = Out4(outs==0); 
%     
%     % split by stimulus type
%     fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
%     for s = 1:length(stmtypes)
%         idx = Lfps_all.stmtype == s;
%         if s > 1
%             continue
%         end
%         % baseline 
%         for f = 1:length(fields)
%             Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
% %             iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
%         end
%         
%         Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
% %         iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
%         
%         % autosave
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_ps/Lfps_mp_' stmtypes{s} '_base.mat'], 'Lfps', '-v7.3')
% %         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_ps/iLfps_mp_' stmtypes{s} '_base.mat'], 'iLfps', '-v7.3')
%         
% %         % drug
% %         for f = 1:length(fields)
% %             Lfps.(fields{f}) = Lfps_dr.(fields{f})(idx);
% %             iLfps.(fields{f}) = iLfps_dr.(fields{f})(idx);
% %         end
% %         
% %         Lfps.LFP_prepro = Lfps_dr.LFP_prepro(idx);
% %         iLfps.LFP_prepro = iLfps_dr.LFP_prepro(idx);
% %         
% %         % autosave
% %         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_ps/Lfps_mp_' stmtypes{s} '_drug.mat'], 'Lfps', '-v7.3')
% %         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_ps/iLfps_mp_' stmtypes{s} '_drug.mat'], 'iLfps', '-v7.3')
% %         
%         disp([stmtypes{s} ': lfps saved!'])
%     end
% end
% 
% % basic LFP analysis =========================
% if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_mp_nothin_ps'))
%     % stimulus types
%     stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
%    
%     % load data and lists
%     loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];    
%     loadpath_mp = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/MPprepro/'];
%     a = load([loadpath 'lfplist.mat'], 'lfplist');  
%     b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
%     c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
%     incl_i = c.incl_i;
%     list_RC = b.list_RC;
%     lfplist = a.lfplist;
%     N = length(lfplist);
%     
%     % initialization
%     Out1 = cell(1, N); 
% %     Out2 = cell(1, N);
% %     Out3 = cell(1, N); 
% %     Out4 = cell(1, N);
%     goodunit = zeros(1, N); is5ht = zeros(1, N);
%     stmtype = zeros(1, N); animal = zeros(1, N); 
%     parfor i = 1:N            
%         % goodunit
%         if isempty(strfind(lfplist{i}{1}, 'xRC'))
%             goodunit(i) = ismember(i, incl_i);
%             continue
%         else
%             goodunit(i) = ismember(i, list_RC);
%             stmtype(i) = 1;
%         end
%         try
%             % load =======================
%             % baseline ---------------------------
%             d0 = load([loadpath lfplist{i}{1}], 'ex');
%             d3 = load([loadpath_mp 'Trials/' lfplist{i}{1}], 'exn');
% %             d5 = load([loadpath_mp 'iTrials/' lfplist{i}{1}], 'exn');
%             
%             % split trials based on pupil
%             [ex, idx] = ex_spliter(d0.ex);
%             exn3_0 = d3.exn;
%             exn3_0.Trials = exn3_0.Trials(idx{1});
%             exn3_2 = d3.exn;
%             exn3_2.Trials = exn3_2.Trials(idx{2});
% %             exn5_0 = d5.exn;
% %             exn5_0.Trials = exn5_0.Trials(idx{1});
% %             exn5_2 = d5.exn;
% %             exn5_2.Trials = exn5_2.Trials(idx{2});
% %             
%             % analysis
%             Out1{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn3_0, exn3_2, 'LFP_prepro', 0);
% %             Out2{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn5_0, exn5_2, 'iLFP_prepro', 0);
%             
%             % drug -----------------------------------
% %             d2 = load([loadpath lfplist{i}{2}], 'ex');
% %             d4 = load([loadpath_mp 'Trials/' lfplist{i}{2}], 'exn');
% %             d6 = load([loadpath_mp 'iTrials/' lfplist{i}{2}], 'exn');
% %             
% %             % split trials based on pupil
% %             [ex, idx] = ex_spliter(d2.ex);
% %             exn4_0 = d4.exn;
% %             exn4_0.Trials = exn4_0.Trials(idx{1});
% %             exn4_2 = d4.exn;
% %             exn4_2.Trials = exn4_2.Trials(idx{2});
% %             exn6_0 = d6.exn;
% %             exn6_0.Trials = exn6_0.Trials(idx{1});
% %             exn6_2 = d6.exn;
% %             exn6_2.Trials = exn6_2.Trials(idx{2});
% %             
% %             % analysis
% %             Out3{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn4_0, exn4_2, 'LFP_prepro', 0);
% %             Out4{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn6_0, exn6_2, 'iLFP_prepro', 0);
% %             
%             % stimlus type
%             if stmtype(i) == 0
%                 stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
%             end
% 
%             % is mango
%             animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');
% 
%             % is 5HT
%             is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));
% 
%             disp(['session ' num2str(i) ' analyzed!'])
%         catch
%             disp(['session ' num2str(i) ' error'])
%         end
%     end
%     % unit info
% %     outs = cellfun('isempty', Out1) | cellfun('isempty', Out2) | ...
% %         cellfun('isempty', Out3) | cellfun('isempty', Out4);
%     outs = cellfun('isempty', Out1);
%     lfplist = lfplist(outs==0);
%     Lfps_all.lfplist = lfplist;
%     Lfps_all.stmtype = stmtype(outs==0);
%     Lfps_all.animal = animal(outs==0);
%     Lfps_all.is5ht = is5ht(outs==0);
%     Lfps_all.goodunit = goodunit(outs==0);
% 
%     % results from analysis    
% %     iLfps_all = Lfps_all;
% %     Lfps_dr = Lfps_all;
% %     iLfps_dr = Lfps_all;
%     Lfps_all.LFP_prepro = Out1(outs==0); 
% %     iLfps_all.LFP_prepro = Out2(outs==0); 
% %     Lfps_dr.LFP_prepro = Out3(outs==0); 
% %     iLfps_dr.LFP_prepro = Out4(outs==0); 
%     
%     % split by stimulus type
%     fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
%     for s = 1:length(stmtypes)
%         idx = Lfps_all.stmtype == s;
%         if s > 1
%             continue
%         end
%         % baseline 
%         for f = 1:length(fields)
%             Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
% %             iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
%         end
%         
%         Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
% %         iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
%         
%         % autosave
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_ps/Lfps_mp_' stmtypes{s} '_base.mat'], 'Lfps', '-v7.3')
% %         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_ps/iLfps_mp_' stmtypes{s} '_base.mat'], 'iLfps', '-v7.3')
%         
%         % drug
% %         for f = 1:length(fields)
% %             Lfps.(fields{f}) = Lfps_dr.(fields{f})(idx);
% %             iLfps.(fields{f}) = iLfps_dr.(fields{f})(idx);
% %         end
% %         
% %         Lfps.LFP_prepro = Lfps_dr.LFP_prepro(idx);
% %         iLfps.LFP_prepro = iLfps_dr.LFP_prepro(idx);
% %         
% %         % autosave
% %         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_ps/Lfps_mp_' stmtypes{s} '_drug.mat'], 'Lfps', '-v7.3')
% %         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_ps/iLfps_mp_' stmtypes{s} '_drug.mat'], 'iLfps', '-v7.3')
% %         
%         disp([stmtypes{s} ': lfps saved!'])
%     end
% end
% 
% 
% % basic LFP analysis =========================
% if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_mp_sc'))
%     % stimulus types
%     stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
%    
%     % load data and lists
%     loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];    
%     loadpath_mp = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/MPprepro/'];
%     a = load([loadpath 'lfplist.mat'], 'lfplist');  
%     b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
%     c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
%     incl_i = c.incl_i;
%     list_RC = b.list_RC;
%     lfplist = a.lfplist;
%     N = length(lfplist);
%     
%     % initialization
%     Out1 = cell(1, N); 
% %     Out2 = cell(1, N);
% %     Out3 = cell(1, N); 
% %     Out4 = cell(1, N);
%     goodunit = zeros(1, N); is5ht = zeros(1, N);
%     stmtype = zeros(1, N); animal = zeros(1, N); 
%     parfor i = 1:N            
%         % goodunit
%         if isempty(strfind(lfplist{i}{1}, 'xRC'))
%             goodunit(i) = ismember(i, incl_i);
%             continue
%         else
%             goodunit(i) = ismember(i, list_RC);
%             stmtype(i) = 1;
%         end
%         try
%             % load =======================
%             % baseline ---------------------------
%             d0 = load([loadpath lfplist{i}{1}], 'ex');
%             d3 = load([loadpath_mp 'Trials/' lfplist{i}{1}], 'exn');
% %             d5 = load([loadpath_mp 'iTrials/' lfplist{i}{1}], 'exn');
%             
%             % split trials based on sc
%             [ex, idx] = ex_spliter(d0.ex, 'sc');
%             exn3_0 = d3.exn;
%             exn3_0.Trials = exn3_0.Trials(idx{1});
%             exn3_2 = d3.exn;
%             exn3_2.Trials = exn3_2.Trials(idx{2});
% %             exn5_0 = d5.exn;
% %             exn5_0.Trials = exn5_0.Trials(idx{1});
% %             exn5_2 = d5.exn;
% %             exn5_2.Trials = exn5_2.Trials(idx{2});
% %             
%             % analysis
%             Out1{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn3_0, exn3_2, 'LFP_prepro');
% %             Out2{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn5_0, exn5_2, 'iLFP_prepro');
%             
%             % drug -----------------------------------
% %             d2 = load([loadpath lfplist{i}{2}], 'ex');
% %             d4 = load([loadpath_mp 'Trials/' lfplist{i}{2}], 'exn');
% %             d6 = load([loadpath_mp 'iTrials/' lfplist{i}{2}], 'exn');
% %             
% %             % split trials based on sc
% %             [ex, idx] = ex_spliter(d2.ex, 'sc');
% %             exn4_0 = d4.exn;
% %             exn4_0.Trials = exn4_0.Trials(idx{1});
% %             exn4_2 = d4.exn;
% %             exn4_2.Trials = exn4_2.Trials(idx{2});
% %             exn6_0 = d6.exn;
% %             exn6_0.Trials = exn6_0.Trials(idx{1});
% %             exn6_2 = d6.exn;
% %             exn6_2.Trials = exn6_2.Trials(idx{2});
% %             
% %             % analysis
% %             Out3{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn4_0, exn4_2, 'LFP_prepro');
% %             Out4{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn6_0, exn6_2, 'iLFP_prepro');
% %             
%             % stimlus type
%             if stmtype(i) == 0
%                 stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
%             end
% 
%             % is mango
%             animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');
% 
%             % is 5HT
%             is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));
% 
%             disp(['session ' num2str(i) ' analyzed!'])
%         catch
%             disp(['session ' num2str(i) ' error'])
%         end
%     end
%     % unit info
% %     outs = cellfun('isempty', Out1) | cellfun('isempty', Out2) | ...
% %         cellfun('isempty', Out3) | cellfun('isempty', Out4);
%     outs = cellfun('isempty', Out1);
%     lfplist = lfplist(outs==0);
%     Lfps_all.lfplist = lfplist;
%     Lfps_all.stmtype = stmtype(outs==0);
%     Lfps_all.animal = animal(outs==0);
%     Lfps_all.is5ht = is5ht(outs==0);
%     Lfps_all.goodunit = goodunit(outs==0);
% 
%     % results from analysis    
% %     iLfps_all = Lfps_all;
% %     Lfps_dr = Lfps_all;
% %     iLfps_dr = Lfps_all;
%     Lfps_all.LFP_prepro = Out1(outs==0); 
% %     iLfps_all.LFP_prepro = Out2(outs==0); 
% %     Lfps_dr.LFP_prepro = Out3(outs==0); 
% %     iLfps_dr.LFP_prepro = Out4(outs==0); 
%     
%     % split by stimulus type
%     fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
%     for s = 1:length(stmtypes)
%         idx = Lfps_all.stmtype == s;
%         if s > 1
%             continue
%         end
%         % baseline 
%         for f = 1:length(fields)
%             Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
% %             iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
%         end
%         
%         Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
% %         iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
%         
%         % autosave
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_sc/Lfps_mp_' stmtypes{s} '_base.mat'], 'Lfps', '-v7.3')
% %         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_sc/iLfps_mp_' stmtypes{s} '_base.mat'], 'iLfps', '-v7.3')
%         
%         % drug
% %         for f = 1:length(fields)
% %             Lfps.(fields{f}) = Lfps_dr.(fields{f})(idx);
% %             iLfps.(fields{f}) = iLfps_dr.(fields{f})(idx);
% %         end
% %         
% %         Lfps.LFP_prepro = Lfps_dr.LFP_prepro(idx);
% %         iLfps.LFP_prepro = iLfps_dr.LFP_prepro(idx);
% %         
% %         % autosave
% %         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_sc/Lfps_mp_' stmtypes{s} '_drug.mat'], 'Lfps', '-v7.3')
% %         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_sc/iLfps_mp_' stmtypes{s} '_drug.mat'], 'iLfps', '-v7.3')
%         
%         disp([stmtypes{s} ': lfps saved!'])
%     end
% end

% basic LFP analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'lfps_pair_mp_nothin_sc'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];    
    loadpath_mp = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/MPprepro_RC/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); 
%     Out2 = cell(1, N);
%     Out3 = cell(1, N); 
%     Out4 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            goodunit(i) = ismember(i, incl_i);
            continue
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            % load =======================
            % baseline ---------------------------
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            d3 = load([loadpath_mp 'Trials/' lfplist{i}{1}], 'exn');
%             d5 = load([loadpath_mp 'iTrials/' lfplist{i}{1}], 'exn');
            
            % split trials based on sc
            [ex0, ex2, sidx1, sidx2] = ex_spliter_es(d0.ex, es, thre, [0.8 0]);
            exn3_0 = d3.exn;
            exn3_0.Trials = exn3_0.Trials(sidx1);
            exn3_2 = d3.exn;
            exn3_2.Trials = exn3_2.Trials(sidx2);
%             exn5_0 = d5.exn;
%             exn5_0.Trials = exn5_0.Trials(idx{1});
%             exn5_2 = d5.exn;
%             exn5_2.Trials = exn5_2.Trials(idx{2});
            
            % analysis
            Out1{i} = pair_stmLFP4mp(ex0, ex2, exn3_0, exn3_2, 'LFP_prepro', 0);
%             Out2{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn5_0, exn5_2, 'iLFP_prepro', 0);
            
            % drug -----------------------------------
%             d2 = load([loadpath lfplist{i}{2}], 'ex');
%             d4 = load([loadpath_mp 'Trials/' lfplist{i}{2}], 'exn');
%             d6 = load([loadpath_mp 'iTrials/' lfplist{i}{2}], 'exn');
%             
%             % split trials based on sc
%             [ex, idx] = ex_spliter(d2.ex, 'sc');
%             exn4_0 = d4.exn;
%             exn4_0.Trials = exn4_0.Trials(idx{1});
%             exn4_2 = d4.exn;
%             exn4_2.Trials = exn4_2.Trials(idx{2});
%             exn6_0 = d6.exn;
%             exn6_0.Trials = exn6_0.Trials(idx{1});
%             exn6_2 = d6.exn;
%             exn6_2.Trials = exn6_2.Trials(idx{2});
%             
%             % analysis
%             Out3{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn4_0, exn4_2, 'LFP_prepro', 0);
%             Out4{i} = pair_stmLFP4mp(ex{1}, ex{2}, exn6_0, exn6_2, 'iLFP_prepro', 0);
%             
            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
%     outs = cellfun('isempty', Out1) | cellfun('isempty', Out2) | ...
%         cellfun('isempty', Out3) | cellfun('isempty', Out4);
    outs = cellfun('isempty', Out1);
    lfplist = lfplist(outs==0);
    Lfps_all.lfplist = lfplist;
    Lfps_all.stmtype = stmtype(outs==0);
    Lfps_all.animal = animal(outs==0);
    Lfps_all.is5ht = is5ht(outs==0);
    Lfps_all.goodunit = goodunit(outs==0);

    % results from analysis    
%     iLfps_all = Lfps_all;
%     Lfps_dr = Lfps_all;
%     iLfps_dr = Lfps_all;
    Lfps_all.LFP_prepro = Out1(outs==0); 
%     iLfps_all.LFP_prepro = Out2(outs==0); 
%     Lfps_dr.LFP_prepro = Out3(outs==0); 
%     iLfps_dr.LFP_prepro = Out4(outs==0); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = Lfps_all.stmtype == s;
        if s > 1
            continue
        end
        % baseline 
        for f = 1:length(fields)
            Lfps.(fields{f}) = Lfps_all.(fields{f})(idx);
%             iLfps.(fields{f}) = iLfps_all.(fields{f})(idx);
        end
        
        Lfps.LFP_prepro = Lfps_all.LFP_prepro(idx);
%         iLfps.LFP_prepro = iLfps_all.LFP_prepro(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_sc/Lfps_mp_' stmtypes{s} '_base.mat'], 'Lfps', '-v7.3')
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_sc/iLfps_mp_' stmtypes{s} '_base.mat'], 'iLfps', '-v7.3')
        
        % drug
%         for f = 1:length(fields)
%             Lfps.(fields{f}) = Lfps_dr.(fields{f})(idx);
%             iLfps.(fields{f}) = iLfps_dr.(fields{f})(idx);
%         end
%         
%         Lfps.LFP_prepro = Lfps_dr.LFP_prepro(idx);
%         iLfps.LFP_prepro = iLfps_dr.LFP_prepro(idx);
%         
%         % autosave
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_sc/Lfps_mp_' stmtypes{s} '_drug.mat'], 'Lfps', '-v7.3')
%         save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Lfps_pair_nothin_sc/iLfps_mp_' stmtypes{s} '_drug.mat'], 'iLfps', '-v7.3')
        
        disp([stmtypes{s} ': lfps saved!'])
    end
end


% IRASA across sessions =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'Irasas'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/IRASAprepro/'];
    a = load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
%     c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
%     incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); Out2 = cell(1, N); 
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            continue
%             goodunit(i) = ismember(i, incl_i);
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            d0 = load([loadpath lfplist{i}{1}], 'exn');
            d2 = load([loadpath lfplist{i}{2}], 'exn');
            
%             % stimlus type
%             if stmtype(i) == 0
%                 if ismember(1, contains(lfplist{i}{1}, '.OR'))
%                     stmtype(i) = 2;
%                 elseif ismember(1, contains(lfplist{i}{1}, '.CO'))
%                     stmtype(i) = 3;
%                 elseif ismember(1, contains(lfplist{i}{1}, '.SF'))
%                     stmtype(i) = 4;
%                 elseif ismember(1, contains(lfplist{i}{1}, '.SZ'))
%                     stmtype(i) = 5;  
%                 end
%             end

            [Out1{i}, Out2{i}] = irasa_pair({d0.exn, d2.exn}, stmtype(i));
            
            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
    Irasas_all.lfplist = lfplist(oks);
    Irasas_all.stmtype = stmtype(oks);
    Irasas_all.animal = animal(oks);
    Irasas_all.is5ht = is5ht(oks);
    Irasas_all.goodunit = goodunit(oks);

    % results from analysis
    Irasas_alli = Irasas_all;
    Irasas_all.pairs = Out1(oks);     
    Irasas_alli.pairs = Out2(oks); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit', 'pairs'};
    for s = 1:length(stmtypes)
        idx = Irasas_all.stmtype == s;
        
        if s > 1
            continue
        end
        
        for f = 1:length(fields)
            Irasas.(fields{f}) = Irasas_all.(fields{f})(idx);
            iIrasas.(fields{f}) = Irasas_alli.(fields{f})(idx);
        end
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Irasas/Irasas_' stmtypes{s} '.mat'], 'Irasas', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Irasas/iIrasas_' stmtypes{s} '.mat'], 'iIrasas', '-v7.3')
        disp([stmtypes{s} ': Irasas saved!'])
    end
end


% IRASA across sessions =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'Irasas_sc'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/IRASAprepro/'];
    a = load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
%     c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
%     incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    es = 0.6592;
    thre = 0.01;
    wnd = [0.8 0];
    
    % initialization
    Out1 = cell(1, N); Out2 = cell(1, N); 
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            continue
%             wnd = [0.2 0];
%             goodunit(i) = ismember(i, incl_i);
        else
%             wnd = [0.8 0];
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            d0 = load([loadpath lfplist{i}{1}], 'exn');
%             d2 = load([loadpath lfplist{i}{2}], 'exn');
            
%             % stimlus type
%             extype = 'or';
%             if stmtype(i) == 0
%                 if ismember(1, contains(lfplist{i}{1}, '.OR'))
%                     stmtype(i) = 2;
%                 elseif ismember(1, contains(lfplist{i}{1}, '.CO'))
%                     stmtype(i) = 3;
%                     extype = 'co';
%                 elseif ismember(1, contains(lfplist{i}{1}, '.SF'))
%                     stmtype(i) = 4;
%                     extype = 'sf';
%                 elseif ismember(1, contains(lfplist{i}{1}, '.SZ'))
%                     stmtype(i) = 5;  
%                     extype = 'sz';
%                 end
%             end

            % split trials based on spike counts
            try
                [ex0, ex2] = ex_spliter_es(d0.exn, es, thre, wnd);
            catch
                % split trials based on spike counts
                ex = ex_spliter(d0.exn, 'sc');
                ex0 = ex{2};
                ex2 = ex{1};
            end
            
            % main analysis
            [Out1{i}, Out2{i}] = irasa_pair({ex0, ex2}, stmtype(i));
            
            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
    Irasas_all.lfplist = lfplist(oks);
    Irasas_all.stmtype = stmtype(oks);
    Irasas_all.animal = animal(oks);
    Irasas_all.is5ht = is5ht(oks);
    Irasas_all.goodunit = goodunit(oks);

    % results from analysis
    Irasas_alli = Irasas_all;
    Irasas_all.pairs = Out1(oks);     
    Irasas_alli.pairs = Out2(oks); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit', 'pairs'};
    for s = 1:length(stmtypes)
        idx = Irasas_all.stmtype == s;
        
        if s > 1
            continue
        end
        
        for f = 1:length(fields)
            Irasas.(fields{f}) = Irasas_all.(fields{f})(idx);
            iIrasas.(fields{f}) = Irasas_alli.(fields{f})(idx);
        end
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Irasas_sc/Irasas_' stmtypes{s} '.mat'], 'Irasas', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Irasas_sc/iIrasas_' stmtypes{s} '.mat'], 'iIrasas', '-v7.3')
        disp([stmtypes{s} ': Irasas saved!'])
    end
end

% IRASA across sessions =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'Irasas_ps'))
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/IRASAprepro/'];
    a = load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
%     c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
%     incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); Out2 = cell(1, N); 
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            continue
%             goodunit(i) = ismember(i, incl_i);
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            d0 = load([loadpath lfplist{i}{1}], 'exn');
%             d2 = load([loadpath lfplist{i}{2}], 'exn');
            
%             % stimlus type
%             extype = 'or';
%             if stmtype(i) == 0
%                 if ismember(1, contains(lfplist{i}{1}, '.OR'))
%                     stmtype(i) = 2;
%                 elseif ismember(1, contains(lfplist{i}{1}, '.CO'))
%                     stmtype(i) = 3;
%                     extype = 'co';
%                 elseif ismember(1, contains(lfplist{i}{1}, '.SF'))
%                     stmtype(i) = 4;
%                     extype = 'sf';
%                 elseif ismember(1, contains(lfplist{i}{1}, '.SZ'))
%                     stmtype(i) = 5;  
%                     extype = 'sz';
%                 end
%             end

            % split trials based on spike counts
            ex = ex_spliter(d0.exn, 'pupil', 2, 'or');
            
            % main analysis
            [Out1{i}, Out2{i}] = irasa_pair({ex{1}, ex{2}}, stmtype(i));
            
            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
    Irasas_all.lfplist = lfplist(oks);
    Irasas_all.stmtype = stmtype(oks);
    Irasas_all.animal = animal(oks);
    Irasas_all.is5ht = is5ht(oks);
    Irasas_all.goodunit = goodunit(oks);

    % results from analysis
    Irasas_alli = Irasas_all;
    Irasas_all.pairs = Out1(oks);     
    Irasas_alli.pairs = Out2(oks); 
    
    % split by stimulus type
    fields = {'lfplist', 'stmtype', 'animal', 'is5ht', 'goodunit', 'pairs'};
    for s = 1:length(stmtypes)
        idx = Irasas_all.stmtype == s;
        
        if s > 1
            continue
        end
        
        for f = 1:length(fields)
            Irasas.(fields{f}) = Irasas_all.(fields{f})(idx);
            iIrasas.(fields{f}) = Irasas_alli.(fields{f})(idx);
        end
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Irasas_ps/Irasas_' stmtypes{s} '.mat'], 'Irasas', '-v7.3')
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/Irasas_ps/iIrasas_' stmtypes{s} '.mat'], 'iIrasas', '-v7.3')
        disp([stmtypes{s} ': Irasas saved!'])
    end
end

% % pairtypes = {'drug', 'sc', 'ps', 'ps_drug'};
% % pairtypes = {'sc', 'ps', 'ps_drug'};
% % pairtypes = {'ps', 'ps_drug'};
% % pairtypes = {'sc'};
% if sum(strcmp(type, 'all')) || sum(strcmp(type,  'basic')
%     % load data and lists
%     a = load([mypath '/Katsuhisa/serotonin_project/dataset/Data/exinfo.mat'], 'exinfo');
%     b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
%     c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
%     exinfo = a.exinfo;
%     list_RC = b.list_RC;
%     incl_i = c.incl_i;
%     % spk and LFP analysis
%     for i = 1:length(pairtypes)
%         for k = 1:2
%             addLFP(exinfo, {list_RC, incl_i}, pairtypes{i}, k-1);
%         end
%     end
% end

% interaction analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'interaction'))
    % load RC data
    load([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/rcdat.mat'], 'dat')
    % interaction between pupil size and 5HT on neural responses
    interaction_lfp_dataset(dat);
end

%  HMM analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'hmm'))    
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); Out2 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            continue
%             goodunit(i) = ismember(i, incl_i);
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            % load =======================
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            d2 = load([loadpath lfplist{i}{2}], 'ex');
            
            % analysis
            Out1{i} = fitHMM_dataset(d0.ex, [0.2, 0.01], 'or_seq', 'oSpikes', [0.8 0]);
            Out2{i} = fitHMM_dataset(d2.ex, [0.2, 0.01], 'or_seq', 'oSpikes', [0.8 0]);
            clc
            
            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
    HMM_all.list = lfplist(oks);
    HMM_all.stmtype = stmtype(oks);
    HMM_all.animal = animal(oks);
    HMM_all.is5ht = is5ht(oks);
    HMM_all.goodunit = goodunit(oks);

    HMM_all.cond(1).results = Out1(oks); 
    HMM_all.cond(2).results = Out2(oks); 
    
    % split by stimulus type
    fields = {'list', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = HMM_all.stmtype == s;
        
        if s > 1
            continue
        end
        
        for f = 1:length(fields)
            hmms.(fields{f}) = HMM_all.(fields{f})(idx);
        end
        
        hmms.cond(1).results = HMM_all.cond(1).results(idx);
        hmms.cond(2).results = HMM_all.cond(2).results(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/hmm/hmms_' stmtypes{s} '.mat'], 'hmms', '-v7.3')
        disp([stmtypes{s} ': hmms saved!'])
    end
end

% Conductance-based neural-encoding model analysis =========================
if sum(strcmp(type, 'all')) || sum(strcmp(type,  'CBEM'))
    addpath(genpath([mypath '/Katsuhisa/code/integrated/CBEM/code']))
    
    % stimulus types
    stmtypes = {'rc', 'or', 'co', 'sf', 'sz'};
   
    % load data and lists
    loadpath = [mypath '/Katsuhisa/serotonin_project/LFP_project/Data/LFPprepro/'];
    a = load([loadpath 'lfplist.mat'], 'lfplist');  
    b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
    c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
    incl_i = c.incl_i;
    list_RC = b.list_RC;
    lfplist = a.lfplist;
    N = length(lfplist);
    
    % initialization
    Out1 = cell(1, N); Out2 = cell(1, N);
    goodunit = zeros(1, N); is5ht = zeros(1, N);
    stmtype = zeros(1, N); animal = zeros(1, N); 
    parfor i = 1:N            
        % goodunit
        if isempty(strfind(lfplist{i}{1}, 'xRC'))
            continue
%             goodunit(i) = ismember(i, incl_i);
        else
            goodunit(i) = ismember(i, list_RC);
            stmtype(i) = 1;
        end
        try
            % load =======================
            d0 = load([loadpath lfplist{i}{1}], 'ex');
            d2 = load([loadpath lfplist{i}{2}], 'ex');
            
            % analysis
            Out1{i} = fitCBEM(d0.ex, 0.001, 'or_seq', 'Spikes', 0);
            Out2{i} = fitCBEM(d2.ex, 0.001, 'or_seq', 'Spikes', 0);
            clc
            
            % stimlus type
            if stmtype(i) == 0
                stmtype(i) = find(strcmp(stmtypes, d0.ex.exp.e1.type));
            end

            % is mango
            animal(i) = strcmp(lfplist{i}{1}(1:2), 'ma');

            % is 5HT
            is5ht(i) = ismember(1, contains(lfplist{i}{2}, '5HT'));

            disp(['session ' num2str(i) ' analyzed!'])
        catch
            disp(['session ' num2str(i) ' error'])
        end
    end
    % unit info
    oks = ~cellfun('isempty', Out1) & ~cellfun('isempty', Out2);
    CBEM_all.list = lfplist(oks);
    CBEM_all.stmtype = stmtype(oks);
    CBEM_all.animal = animal(oks);
    CBEM_all.is5ht = is5ht(oks);
    CBEM_all.goodunit = goodunit(oks);

    CBEM_all.cond(1).results = Out1(oks); 
    CBEM_all.cond(2).results = Out2(oks); 
    
    % split by stimulus type
    fields = {'list', 'stmtype', 'animal', 'is5ht', 'goodunit'};
    for s = 1:length(stmtypes)
        idx = CBEM_all.stmtype == s;
        
        if s > 1
            continue
        end
        
        for f = 1:length(fields)
            cbem.(fields{f}) = CBEM_all.(fields{f})(idx);
        end
        
        cbem.cond(1).results = CBEM_all.cond(1).results(idx);
        cbem.cond(2).results = CBEM_all.cond(2).results(idx);
        
        % autosave
        save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/CBEM/cbem_' stmtypes{s} '.mat'], 'cbem', '-v7.3')
        disp([stmtypes{s} ': cbem saved!'])
    end
end