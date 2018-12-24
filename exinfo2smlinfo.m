function smlinfo = exinfo2smlinfo(exinfo)
% extract relevant info from exinfo

% path =======================================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end

% initialize
smlinfo.columnames = {'ismango', 'is5ht', 'isgood', 'stmtype', 'additive change', 'gain change', ...
    'wavewidth', 'RFx', 'RFy', 'noise corr base', 'noise corr drug', 'fano factor base', 'fano factor drug', ...
    'relative rate'};
smlinfo.paramat = nan(length(exinfo), length(smlinfo.columnames));
stmtypes = {'or', 'co', 'sz', 'sf'};
b = load([mypath '/Corinna/SharedCode/Katsu/list_RC.mat'], 'list_RC'); 
c = load([mypath '/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat'], 'incl_i'); 
incl_i = c.incl_i;
list_RC = b.list_RC;
wnd = [0.2 0];

% extract info
for i = 1:length(exinfo)
  % parameter matrix =========================
  % is mango
  smlinfo.paramat(i, 1) = 1*strcmp(exinfo(i).monkey, 'ma');
  
  % is5ht
  smlinfo.paramat(i, 2) = 1*strcmp(exinfo(i).drugname, '5HT');
  
  % is goodunit
  if exinfo(i).isRC
     smlinfo.paramat(i, 3) = ismember(i, list_RC);
  else
     smlinfo.paramat(i, 3) = ismember(i, incl_i);
  end
  
  % stimulus type
  if exinfo(i).isRC
    smlinfo.paramat(i, 4) = 0;
  else
    smlinfo.paramat(i, 4) = find(contains(stmtypes, exinfo(i).param1));
  end
  
  % additive change
  smlinfo.paramat(i, 5) = exinfo(i).yoff;
  
  % gain change
  smlinfo.paramat(i, 6) = exinfo(i).gslope;
  
  % wavewidth
  smlinfo.paramat(i, 7) = mean(exinfo(i).wdt);
  
  % RF x
  smlinfo.paramat(i, 8) = exinfo(i).RFwx(1);
  
  % RF y
  smlinfo.paramat(i, 9) = exinfo(i).RFwy(1);
  
  % noise correlation & fano factor
  if isempty(exinfo(i).rsc) || isempty(exinfo(i).rsc_drug)
      fname = strrep(exinfo(i).fname, '\', '/');
      fname = [mypath fname(3:end)];
      [smlinfo.paramat(i, 10), smlinfo.paramat(i, 12)] = getNC(fname, wnd);
      fname = strrep(exinfo(i).fname_drug, '\', '/');
      fname = [mypath fname(3:end)];
      [smlinfo.paramat(i, 11), smlinfo.paramat(i, 13)] = getNC(fname, wnd);
  else
      smlinfo.paramat(i, 10) = exinfo(i).rsc;
      smlinfo.paramat(i, 11) = exinfo(i).rsc_drug;
      smlinfo.paramat(i, 12) = nanmean(exinfo(i).ff.classic.ff);
      smlinfo.paramat(i, 13) = nanmean(exinfo(i).ff_drug.classic.ff);
  end
    
  % relative rate
  if isempty(exinfo(i).nonparam_ratio)
     smlinfo.paramat(i, 14) = mean(exinfo(i).fitparam_drug.val.mn)/mean(exinfo(i).fitparam.val.mn);
  else
     smlinfo.paramat(i, 14) = exinfo(i).nonparam_ratio;
  end
  
  % others
  smlinfo.fitparam{i} = exinfo(i).fitparam;
  smlinfo.fitparam_drug{i} = exinfo(i).fitparam_drug;
  smlinfo.fname{i} = exinfo(i).fname;
  smlinfo.fname_drug{i} = exinfo(i).fname_drug;
  
  disp(['session ' num2str(i) ' stored!'])
end

% autosave
save([mypath '/Katsuhisa/serotonin_project/LFP_project/Data/smlinfo.mat'], 'smlinfo', '-v7.3')
disp('smlinfo saved!')