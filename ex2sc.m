function sc = ex2sc(ex)
% spike counts
b = 10;
len_tr = length(ex.Trials);
ncol = round((1000/b)*ex.fix.duration);
sc = zeros(len_tr, ncol);
for i = 1:len_tr
    begin = 0;
    for c = 1:ncol
        sc(i,c) = sum(ex.Trials(i).oSpikes > begin & ex.Trials(i).oSpikes <= begin + 0.001*b);
        begin = begin + 0.001*b;
    end
end