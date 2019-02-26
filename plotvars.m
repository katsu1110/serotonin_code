function plotvars(anaT, varname1, varname2)
%%
% scatter for a pair of variables
%

% addpath(genpath('Z:\Katsuhisa\code\integrated\gramm'))
addpath(genpath('Z:\Katsuhisa\code\integrated\matlab_usefulfunc'))

% vals
X = anaT.table(:, contains(anaT.varnames, varname1));
Y = anaT.table(:, contains(anaT.varnames, varname2));

mat = [X, Y, anaT.lists(:,2), anaT.lists(:,3)];
close all;
ndrug = 1;
drugnames = {'FR'};
if strcmp(anaT.pairnames{1,1}, 'base')
    ndrug = 2;
    drugnames = {'NaCl', '5HT'};
else
    mat(:,3) = 0;
end
for d = 1:ndrug
    subplot(1,ndrug,d)
    if strcmp(varname1(1:end-4), varname2(1:end-4))
        unity_scatter(mat(mat(:,3)==d-1, 1), mat(mat(:,3)==d-1, 2), ...
            mat(mat(:,3)==d-1, 4))
    else
        ls_scatter(mat(mat(:,3)==d-1, 1), mat(mat(:,3)==d-1, 2), ...
            mat(mat(:,3)==d-1, 4))
    end
    title(drugnames{d})
    if length(drugnames) > 1
        xlabel(varname1)
        ylabel(varname2)
    else
        xlabel([varname1(1:end-4) ' (high FR)'])
        ylabel([varname2(1:end-4) ' (low FR)'])
    end
end
set(gcf, 'position', [282   340   560   225])
set(gcf, 'Name', [varname1 ' vs ' varname2], 'NumberTitle', 'off')

% stats    
if ndrug > 1
    kaki_nacl = mat(mat(:,3)==0 & mat(:,4)==0, 1:2);
    kaki_5ht = mat(mat(:,3)==1 & mat(:,4)==0, 1:2);
    mango_nacl = mat(mat(:,3)==0 & mat(:,4)==1, 1:2);
    mango_5ht = mat(mat(:,3)==1 & mat(:,4)==1, 1:2);
    stats_kaki = pair_tests(kaki_nacl, kaki_5ht);
    disp(['Kaki (NaCl; n=' num2str(size(kaki_nacl, 1)) ') ------------------------'])
    disp(stats_kaki.pair(1).table)
    disp(['Kaki (5HT; n=' num2str(size(kaki_5ht, 1)) ') ------------------------'])
    disp(stats_kaki.pair(2).table)
    disp(['Kaki (NaCl vs 5HT; n=' num2str(size(kaki_5ht, 1)+size(kaki_nacl, 1)) ') ---------'])
    disp(stats_kaki.table)
    stats_mango = pair_tests(mango_nacl, mango_5ht);
    disp(['Mango (NaCl; n=' num2str(size(mango_nacl, 1)) ') ------------------------'])
    disp(stats_mango.pair(1).table)
    disp(['Mango (5HT; n=' num2str(size(mango_5ht, 1)) ') ------------------------'])
    disp(stats_mango.pair(2).table)
    disp(['Mango (NaCl vs 5HT; n=' num2str(size(mango_5ht, 1)+size(mango_nacl, 1)) ') ---------'])
    disp(stats_mango.table)
    stats_all = pair_tests([kaki_nacl; mango_nacl], [kaki_5ht; mango_5ht]);
    disp(['Both (NaCl; n=' num2str(size(mango_nacl, 1)+size(kaki_nacl, 1)) ') ------------------------'])
    disp(stats_all.pair(1).table)
    disp(['Both (5HT; n=' num2str(size(mango_5ht, 1)+size(kaki_5ht, 1)) ') ------------------------'])
    disp(stats_all.pair(2).table)
    disp(['Both (NaCl vs 5HT; n=' num2str(size(mango_nacl, 1)+size(mango_5ht, 1)+size(kaki_nacl, 1)+size(kaki_5ht, 1)) ') ---------'])
    disp(stats_all.table)
end

% % visualize
% close all;
% g = gramm('x', X, 'y', Y, 'color', anaT.lists(:,3));
% g.facet_grid([], anaT.lists(:,2));
% g.geom_point();
% g.stat_glm();
% g.set_names('column','drug','x',varname1,'y',varname2,'color','animals');
% g.set_title('correlation');
% figure('Position',[100 100 800 400]);
% 
% g.draw();
% 
% % stats    
% kaki_nacl = [g.results.geom_point_handle(1).XData', ...
%    g.results.geom_point_handle(1).YData'];
% kaki_5ht = [g.results.geom_point_handle(3).XData', ...
%    g.results.geom_point_handle(3).YData'];
% mango_nacl = [g.results.geom_point_handle(2).XData', ...
%    g.results.geom_point_handle(2).YData'];
% mango_5ht = [g.results.geom_point_handle(4).XData', ...
%    g.results.geom_point_handle(4).YData'];
% stats_kaki = pair_tests(kaki_nacl, kaki_5ht);
% disp(['Kaki (NaCl; n=' num2str(size(kaki_nacl, 1)) ') ------------------------'])
% disp(stats_kaki.pair(1).table)
% disp(['Kaki (5HT; n=' num2str(size(kaki_5ht, 1)) ') ------------------------'])
% disp(stats_kaki.pair(2).table)
% disp(['Kaki (NaCl vs 5HT; n=' num2str(size(kaki_5ht, 1)+size(kaki_nacl, 1)) ') ---------'])
% disp(stats_kaki.table)
% stats_mango = pair_tests(mango_nacl, mango_5ht);
% disp(['Mango (NaCl; n=' num2str(size(mango_nacl, 1)) ') ------------------------'])
% disp(stats_mango.pair(1).table)
% disp(['Mango (5HT; n=' num2str(size(mango_5ht, 1)) ') ------------------------'])
% disp(stats_mango.pair(2).table)
% disp(['Mango (NaCl vs 5HT; n=' num2str(size(mango_5ht, 1)+size(mango_nacl, 1)) ') ---------'])
% disp(stats_mango.table)
% stats_all = pair_tests([kaki_nacl; mango_nacl], [kaki_5ht; mango_5ht]);
% disp(['Both (NaCl; n=' num2str(size(mango_nacl, 1)+size(kaki_nacl, 1)) ') ------------------------'])
% disp(stats_all.pair(1).table)
% disp(['Both (5HT; n=' num2str(size(mango_5ht, 1)+size(kaki_5ht, 1)) ') ------------------------'])
% disp(stats_all.pair(2).table)
% disp(['Both (NaCl vs 5HT; n=' num2str(size(mango_nacl, 1)+size(mango_5ht, 1)+size(kaki_nacl, 1)+size(kaki_5ht, 1)) ') ---------'])
% disp(stats_all.table)