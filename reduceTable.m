function anaTred = reduceTable(anaT, focus)
%%
% reduce the dimension of the table
%

% sets ====================
if nargin < 2; focus = 'clustering'; end
switch focus
    case 'clustering'
%         vars = {'r rate', 'd noise corr', 'd selectivity', 'd snr2', ...
%             'additive change', 'gain change', 'wavewidth', 'RFx', 'RFy', 'd ff'};
%         vars = {'additive change', 'gain change', 'wavewidth', 'RFx', 'RFy'};
        vars = {'additive change', 'gain change', 'wavewidth', 'RFx', 'RFy', 'r rate'};
    case 'LFPpow' 
        vars = {'d theta pow', 'd alpha pow', 'd beta pow', 'd gamma pow',...
            'd low-freq pow', 'd broadband pow'};
    case 'LFPres'
        vars = {'d theta res', 'd alpha res', 'd beta res', 'd gamma res',...
            'd low-freq res', 'd broadband res'};
end

ind = contains(anaT.varnames, vars);
anaTred = anaT;
anaTred.table = anaT.table(:, ind);
anaTred.varnames = anaT.varnames(ind);