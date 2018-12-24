function exs = unitpair_preprocess(exs, period, internal, names)
% shared preprocessing across LFP data (after each 'freqAnalysis.m')

%%
% input
if nargin < 2; period = {[0.25, 0.45]}; end
if nargin < 3; internal = 0; end
if nargin < 4
    ex = exs{1};
    if isfield(ex,'fileID')
        filename = ex.fileName;
    elseif isfield(ex,'Header')
        if isfield(ex.Header,'onlineFileName')
            filename = ex.Header.onlineFileName;
        elseif isfield(ex.Header,'Headers')
            filename = ex.Header.fileName;
        else
         filename = ex.Header.fileName;
        end
    else
        filename = ex.fileName;
    end
    names = {filename};
end

if internal == 0
    field = 'LFP_prepro';
elseif internal == 1
    field = 'iLFP_prepro';
end    
params = define_params;

% LFP vector
lenx = length(exs);
v=  [];
for i = 1:lenx
    exs{i}.period = period;
    for k = 1:length(exs{i}.Trials)
        v = [v, exs{i}.Trials(k).(field)(...
            exs{i}.Trials(k).LFP_prepro_time > 0 ...
            & exs{i}.Trials(k).LFP_prepro_time < period{end}(2))];
    end
end

%%
% preprocessing
Fs = 1000; % sampling frequency
folderName = '../../../LFP_project/Data/MPprepro/';
max_iter = 500;
wrap = 1; % if 1, it does not work...
L = 4096; % 2^N
if (period{end}(2) + 0.2)*Fs < 1024
    L = 1024;
elseif (period{end}(2) + 0.2)*Fs < 2048
    L = 2048;
end
signalRange = [1 L];
atomList = [];
me = mean(v(:));
sd = std(v(:));
for i = 1:lenx
    % z-scoring =====================
    ntr = length(exs{i}.Trials);
    exs{i}.FREQ = 0:Fs/L:100;
    lenf = length(exs{i}.FREQ);
    inputSignal = nan(L, ntr, 1);
    idxs = ones(2, ntr);
    for k = 1:ntr
        % z-scoring
        exs{i}.Trials(k).LFP_z = (exs{i}.Trials(k).(field) - me)/sd;
        idxs(1, k) = length(exs{i}.Trials(k).LFP_z);
        
        % padding and mean subtraction for signal processing
        [z, idxs(2, k)] = padding(exs{i}.Trials(k).LFP_z, L);
        inputSignal(:, k, 1) = z - mean(z);
    end
    
    % MP =========================
    tag = names{i};
    
    % format for preprocessing
    importData(inputSignal, folderName, tag, signalRange, Fs);
    
    % perform Gabor decomposition
    runGabor(folderName, tag, L, max_iter);
    
    % signal reconstruction
    minwnd = min(idxs(1, :));
    gaborInfo{1} = getGaborData(folderName, tag, 1);
    for k = 1:ntr
        try
            % signal
            mp_signal = reconstructSignalFromAtomsMPP(...
                gaborInfo{1}{k}.gaborData, L, wrap, atomList);
            mp_signal = mp_signal(idxs(2, k):idxs(2,k)+minwnd-1);
            if isempty(mp_signal)
                error('')
            end
        catch
            disp(['trial ' num2str(k) ': MP signal reconstruction error. Using LFP_z...'])
            mp_signal = exs{i}.Trials(k).LFP_z(1:minwnd);
        end
        exs{i}.Trials(k).MPsignal = mp_signal;
        
        % energy
        try
            rEnergy = reconstructEnergyFromAtomsMPP(...
                gaborInfo{1}{k}.gaborData, L, wrap, atomList);
            rEnergy = rEnergy(1:lenf, idxs(2, k):idxs(2,k)+minwnd-1);
        catch
            disp(['trial ' num2str(k) ': MP signal reconstruction error. Using MTM instead...'])
            rEnergy = mtspecgramc(exs{i}.Trials(k).MPsignal, [0.1 0.01], params);
            rEnergy = imresize(rEnergy, [lenf, minwnd]);
        end
        exs{i}.Trials(k).MPenergy = rEnergy;
    end
end


%         % period: baseline, stimulus evoked, sustained
%         for u = 1:3
%             % LFP trace
%             exs{i}.Trials(k).period(u).LFP_z_time = time(time>=exs{i}.period{u}(1) & time <= exs{i}.period{u}(2));
%             exs{i}.Trials(k).period(u).LFP_z = exs{i}.Trials(k).LFP_z(time>=exs{i}.period{u}(1) & time <= exs{i}.period{u}(2));
%         end
%     end
% end

function [b, pre] = padding(a, L)
% add padding for signal processing
lena = length(a);
b = zeros(1, L);
pre = floor((L - lena)/2) + 1;
b(1:pre-1) = a(1);
b(pre:pre+lena-1) = a;
b(pre+lena:end) = a(end);
