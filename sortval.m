function idx = sortval(val, nsplit)
[~, sortedidx] = sort(val);
d = round(length(val)/nsplit);
idx = cell(1, nsplit);
begin = 1;
for i = 1:nsplit-1
    idx{i} = sortedidx(begin:begin + d - 1);
    begin = begin + d;
end
idx{nsplit} = sortedidx(begin:end);