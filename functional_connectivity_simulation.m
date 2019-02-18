function functional_connectivity_simulation(repeat)
%%
% what would be a good way to use a supervising learning to compare two
% conditions in terms of functional connectivity?
%
if nargin < 1; repeat = 1000; end
rng default  % For reproducibility
metrics = zeros(repeat, 3);   
for r = 1:repeat
    sd = rand(1,2);
    sd = sort(sd, 'descend');
    
    % signals
    mu = [2 3];
    sigma = [1 sd(1); sd(1) 1];
    base = mvnrnd(mu,sigma,1000);
    sigma = [1 sd(2); sd(2) 1];
    drug = mvnrnd(mu,sigma,1000);

    % check correlation
    r1 = corrcoef(base(:,1), base(:,2));
    r2 = corrcoef(drug(:,1), drug(:,2));
    metrics(r, 3) = r1(1, 2) - r2(1, 2);
    disp(['repeat ' num2str(r) '; correlation in base: ' num2str(r1(1,2))])
    disp(['repeat ' num2str(r) '; correlation in drug: ' num2str(r2(1,2))])

    % compute cross-validation score
    cv = 10;
    indices = crossvalind('Kfold', size(base, 1), cv);
    for i = 1:cv
        test = (indices == i); 
        train = ~test;

        % model fitting --- separately ---
        beta_base = glmfit(base(train, 1), base(train, 2));
        ypred = glmval(beta_base, base(test, 1), 'identity');
        r_base = corrcoef(ypred, base(test, 2));
        beta_drug = glmfit(drug(train, 1), drug(train, 2));
        ypred = glmval(beta_drug, drug(test, 1), 'identity');
        r_drug = corrcoef(ypred, drug(test, 2));
        metrics(r, 1) = metrics(r, 1) + (r_base(1,2) - r_drug(1,2));

        % model fitting --- train on base, test on drug ---
        ypred = glmval(beta_base, drug(test, 1), 'identity');
        r_drug = corrcoef(ypred, drug(test, 2));
        metrics(r, 2) = metrics(r, 2) + (r_base(1,2) - r_drug(1,2));
    end
end
metrics(:, 1:2) = metrics(:, 1:2)./cv;

close all;
figure;
plot(metrics(:, 1), metrics(:, 3), 'o', 'color', 'c')
hold on;
plot(metrics(:, 2), metrics(:, 3), 'o', 'color', 'm')
set(gca, 'box', 'off', 'tickdir', 'out')

