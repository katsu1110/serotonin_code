function Lfps = co_control
%
% automatic analysis using co control experiments
%

% unit numbers
cos = {'0367', '0368', '0369', '0371', '0372', '0373'};

% cos = {'0368', '0369', '0371', '0372', '0373'};
% cos = {'0371', '0372', '0373'};
%  cos = {'0367', '0368', '0369'};
lenc = length(cos);

% path ==========================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end

addpath(genpath([mypath '/Katsuhisa/serotonin_project/']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

folderpath = '/data/kaki/';
% ext = {'c0', 'c1', 'lfp'};
es = 0.6527; % effect size of the main experiment
Lfps = [];
for c = 1:lenc
    % folder contents
    listings = dir([mypath folderpath cos{c}]);
    lenl = length(listings);
    
    % contcatenation
    conc = 0;
    for l = 1:lenl
        if ismember(1, contains(listings(l).name, '_c1_'))
            if conc==0
                ex = load([mypath folderpath cos{c} '/' listings(l).name], 'ex');
                ex = co_insert(ex.ex);
                ex_c1_all = ex;
                ex = load([mypath folderpath cos{c} '/' strrep(listings(l).name, 'c1', 'c0')], 'ex');
                ex = co_insert(ex.ex);
                ex_c0_all = ex;
                ex = load([mypath folderpath cos{c} '/' strrep(listings(l).name, 'c1', 'lfp')], 'ex');
                ex = co_insert(ex.ex);
                ex_lfp_all = ex;
            else
                ex = load([mypath folderpath cos{c} '/' listings(l).name], 'ex');
                ex = co_insert(ex.ex);
                old_lentr = length(ex_c1_all.Trials);
                new_lentr = length(ex.Trials);
                if isfield(ex.Trials, 'pa')
                    ex.Trials = rmfield(ex.Trials, 'pa');
                end
                ex_c1_all.Trials(old_lentr+1:old_lentr+new_lentr) = ex.Trials;
                ex = load([mypath folderpath cos{c} '/' strrep(listings(l).name, 'c1', 'c0')], 'ex');
                ex = co_insert(ex.ex);
                if isfield(ex.Trials, 'pa')
                    ex.Trials = rmfield(ex.Trials, 'pa');
                end
                ex_c0_all.Trials(old_lentr+1:old_lentr+new_lentr) = ex.Trials; 
                ex = load([mypath folderpath cos{c} '/' strrep(listings(l).name, 'c1', 'lfp')], 'ex');
                ex = co_insert(ex.ex);
                if isfield(ex.Trials, 'pa')
                    ex.Trials = rmfield(ex.Trials, 'pa');
                end
                ex_lfp_all.Trials(old_lentr+1:old_lentr+new_lentr) = ex.Trials;                
            end
            conc = conc + 1;
        end
    end
    
    % autosave 
    prefix = listings(5).name(1:8);
    sorter = listings(5).name(12:19);
    postfix = listings(5).name(25:end);
    ex = ex_c1_all;
    save([mypath folderpath cos{c} '/' prefix 'c1' sorter 'all' postfix], 'ex')
    ex = ex_c0_all;
    save([mypath folderpath cos{c} '/' prefix 'c0' sorter 'all' postfix], 'ex')
    ex = ex_lfp_all;
    save([mypath folderpath cos{c} '/' prefix 'lfp' sorter 'all' postfix], 'ex')

    % LFP analysis
    ex = loadCluster([mypath folderpath cos{c} '/' prefix 'c1' sorter 'all' postfix], 'loadlfp', 1);
    [ex0, ex2] = split_co(ex, es);
    ex2.exp.e1.vals = ex2.Trials(1).co;
    for j = 1:length(ex2.Trials)
        ex2.Trials(j).co = ex0.Trials(1).co;
    end
    lfps = pair_stmLFP(ex0, ex2, 'LFP_prepro', 0);
    Lfps.lfplist{c} = [mypath folderpath cos{c} '/' prefix 'c1' sorter 'all' postfix];
    Lfps.stmtype(c) = 3;
    Lfps.animal(c) = 0;
    Lfps.is5ht(c) = 0;
    Lfps.goodunit(c) = 1;
    Lfps.LFP_prepro{c} = lfps;
%     lfps.cond(1).si = ex0.si;
%     lfps.cond(2).si = ex2.si;
    save([mypath folderpath cos{c} '/' prefix 'lfps.mat'], 'lfps')
    
    disp([cos{c} ' analyzed!'])
end
save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/control/Lfps_co.mat'], 'Lfps')
disp([cos{c} ' analyzed!'])

function ex = co_insert(ex)
ex.exp.e2 = ex.exp.e1;
ex.exp.e1.type = 'co';
ex.exp.e1.vals = ex.stim.vals.co_range;
for i = 1:length(ex.Trials)
    ex.Trials(i).co = ex.stim.vals.co_range;
end