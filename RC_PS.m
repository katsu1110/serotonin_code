function rc_results = RC_PS(ex0, ex2, ex0_sps, ex0_lps, ...
    ex2_sps, ex2_lps, fields, drugname)
% perform the subspace RC analysis with respect to 5HT and pupil size

rc_results.fields = fields;
% try
    for f = 1:length(fields)
        % replicate Corinna's findings (the effect of 5HT)
        [stmMat, actMat] = ex4RCsub(ex0, 'or', fields{f});
        rc_results.rcsub_base(f).results = reverse_corr_subspace(stmMat, actMat, 300, 100, 0);
        [stmMat, actMat] = ex4RCsub(ex2, 'or', fields{f});
        rc_results.rcsub_drug(f).results = reverse_corr_subspace(stmMat, actMat, 300, 100, 0);
        xv = [rc_results.rcsub_base(f).results.stm.totalcounts];
        yv = [rc_results.rcsub_drug(f).results.stm.totalcounts];
        m = max(xv);
        rc_results.type2reg(f).drug = gmregress(xv/m, yv/m);
        for i = 1:4
            switch i
                case 1
                    exd = ex2_sps;
                    labd = ['pupil small, ' drugname];
                case 2
                    exd = ex2_lps;
                    labd = ['pupil large, ' drugname];
                case 3
                    exd = ex0_sps;
                    labd = 'pupil small, baseline';
                case 4
                    exd = ex0_lps;
                    labd = 'pupil large, baseline';            
            end               
            [stmMat, actMat] = ex4RCsub(exd, 'or', fields{f});
            rc_results.rcsub(f).each(i).results = reverse_corr_subspace(stmMat, actMat, 300, 100, 0);
            rc_results.rcsub(f).each(i).label = labd;
            rc_results.rcsub(f).inter_table_lat(i) = rc_results.rcsub(f).each(i).results.latency;
        end
        % the effect of PS
        ex_sps = concatenate_ex(ex0_sps, ex2_sps);
        [stmMat, actMat] = ex4RCsub(ex_sps, 'or', fields{f});
        rc_results.rcsub_sps(f).results = reverse_corr_subspace(stmMat, actMat, 300, 100, 0);
        ex_lps = concatenate_ex(ex0_lps, ex2_lps);
        [stmMat, actMat] = ex4RCsub(ex_lps, 'or', fields{f});
        rc_results.rcsub_lps(f).results = reverse_corr_subspace(stmMat, actMat, 300, 100, 0);
        xv = [rc_results.rcsub_lps(f).results.stm.totalcounts];
        yv = [rc_results.rcsub_sps(f).results.stm.totalcounts];
        m = max(xv);
        rc_results.type2reg(f).ps = gmregress(xv/m, yv/m);
        % gain or additive change 
        xv = [rc_results.rcsub(f).each(3).results.stm.totalcounts];
        yv = [rc_results.rcsub(f).each(1).results.stm.totalcounts];
        m = max(xv);
        rc_results.type2reg(f).sps_drug = gmregress(xv/m, yv/m);
        xv = [rc_results.rcsub(f).each(4).results.stm.totalcounts];
        yv = [rc_results.rcsub(f).each(2).results.stm.totalcounts];
        m = max(xv);
        rc_results.type2reg(f).lps_drug = gmregress(xv/m, yv/m);    
    end
% catch
%     disp('Provide ex files with the RC experiment.')
%     rc_results = [];
% end

function ex = concatenate_ex(ex0, ex1)
ex = ex0;
len0 = length(ex0.Trials);
len1 = length(ex1.Trials);
ex.Trials(len0+1:len0+len1) = ex1.Trials;