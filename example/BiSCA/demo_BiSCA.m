% demo_ieeg_bixialpha

close all

% Resolve repo root (robust to current folder)
this_file = mfilename('fullpath');
candidate_roots = {
    fileparts(fileparts(fileparts(this_file)));  % 3-deep
    fileparts(fileparts(this_file));             % 2-deep
};

repo_root = [];
for i = 1:numel(candidate_roots)
    if exist(fullfile(candidate_roots{i}, 'setup.m'), 'file') == 2
        repo_root = candidate_roots{i};
        break;
    end
end
assert(~isempty(repo_root), 'Cannot find repository root containing setup.m');

if ~exist('iEEG', 'var')
    filepath = fullfile(repo_root, 'data', 'iEEG', 'iEEG1772FemTposHmInvScMleLog.mat');
    iEEG = load(filepath);
    iEEG = remove_useless_ieeg_fields(iEEG);
end

%%
[subj, ch] = deal(32, 3);
taskname = 'fig1';

[iEEG] = batch_iEEG_hos_bixialpha(iEEG, [subj, ch], repo_root, taskname, 1);

bixialpha = iEEG.subjinfo(subj).bixialpha(ch);
figure(101)
bixialpha.show

%%
function [iEEG] = batch_iEEG_hos_bixialpha(iEEG, idx_subj_ch, repo_root, taskname, ids)
    Fs = 200;

    for id = ids
        isubj = idx_subj_ch(id, 1);
        ich = idx_subj_ch(id, 2);
        disp(['[taskmng]:start task ', num2str(isubj), '_', num2str(ich)]);
        get_freemem(true);
        savepath = fullfile(repo_root, 'result', 'iEEG', 'HoXiAlpha/', ...
            ['W', filesep, taskname, '_cpct', filesep], ...
            ['iEEG_', num2str(isubj), '_', num2str(ich), '.mat']);

        if exist(savepath, "file")

            try
                load(savepath, 'hos', 'nlng', 'bixialpha');
                hos.compact_result();
                nlng.compact_result();
                bixialpha.compact_result();
                iEEG.subjinfo(isubj).bixialpha(ich, :) = bixialpha;

            catch
                warning(['file crash ', num2str(isubj), '_', num2str(ich), ', try to recalculate'])
                iEEG = calc_iEEG_hos_bixialpha(iEEG, repo_root, taskname, Fs, isubj, ich);
            end

        else
            warning(['lost file ', num2str(isubj), '_', num2str(ich), ', try to recalculate'])
            iEEG = calc_iEEG_hos_bixialpha(iEEG, repo_root, taskname, Fs, isubj, ich);
        end

    end

end

%%
function iEEG = calc_iEEG_hos_bixialpha(iEEG, repo_root, taskname, Fs, isubj, ich)

    tic
    x = double(iEEG.subjinfo(isubj).Data_W(:, ich));
    x = x(x ~= 0);
    iEEG.subjinfo(isubj).hos(ich, :) = HOS(x, Fs);
    iEEG.subjinfo(isubj).hos(ich, :).window = 300;
    iEEG.subjinfo(isubj).hos(ich, :).nfft = 300;
    iEEG.subjinfo(isubj).hos(ich, :).frange = [0, 50];

    iEEG.subjinfo(isubj).hos(ich, :).overlap = 0.75;
    iEEG.subjinfo(isubj).hos(ich, :).normalization = 'haubrich'; % haubrich skewness
    % NW=1.5 selected for this dataset
    iEEG.subjinfo(isubj).hos(ich, :).NW = 1.5; %1.5
    iEEG.subjinfo(isubj).hos(ich, :).method = 'pmtm'; % fft pmtm
    iEEG.subjinfo(isubj).hos(ich, :).tf_return_res = false; %true false
    iEEG.subjinfo(isubj).hos(ich, :).tf_return_seg = false; %true false
    iEEG.subjinfo(isubj).hos(ich, :).log_var_s = true;
    iEEG.subjinfo(isubj).hos(ich, :).runall();

    %iEEG info
    iEEG.subjinfo(isubj).hos(ich, :).name = iEEG.subjinfo(isubj).ChName{ich};
    iEEG.subjinfo(isubj).hos(ich, :).info.isubj = isubj;
    iEEG.subjinfo(isubj).hos(ich, :).info.ich = ich;
    iEEG.subjinfo(isubj).hos(ich, :).info.stg = 'W';
    iEEG.subjinfo(isubj).hos(ich, :).info.ptid = iEEG.subjinfo(isubj).PatientId;
    iEEG.subjinfo(isubj).hos(ich, :).info.cname = iEEG.subjinfo(isubj).ChName{ich};
    iEEG.subjinfo(isubj).hos(ich, :).info.ctype = iEEG.subjinfo(isubj).ChannelType{ich};
    iEEG.subjinfo(isubj).hos(ich, :).info.hemis = iEEG.subjinfo(isubj).Hemisphere{ich};
    iEEG.subjinfo(isubj).hos(ich, :).info.cpos = iEEG.subjinfo(isubj).ChannelPositionMid(ich, :);
    iEEG.subjinfo(isubj).hos(ich, :).info.gender = iEEG.subjinfo(isubj).Gender{1};
    iEEG.subjinfo(isubj).hos(ich, :).info.age = iEEG.subjinfo(isubj).AgeAtTimeOfStudy;

    ksall = 0;
    khall = 3;

    if isfield(iEEG.subjinfo(isubj), 'bixialpha') && isnumeric(iEEG.subjinfo(isubj).bixialpha)
        iEEG.subjinfo(isubj).bixialpha = BiXiAlphaBootstrap();
    end

    iEEG.subjinfo(isubj).bixialpha(ich, :) = BiXiAlphaBootstrap(iEEG.subjinfo(isubj).hos(ich, :), ksall, khall);
    iEEG.subjinfo(isubj).bixialpha(ich, :).verbose = false; % true false
    iEEG.subjinfo(isubj).bixialpha(ich, :).multiple_seg = iEEG.subjinfo(isubj).hos(ich, :).tf_return_seg;
    iEEG.subjinfo(isubj).bixialpha(ich, :).runfit;
    iEEG.subjinfo(isubj).bixialpha(ich, :).logging.log = [];
    iEEG.subjinfo(isubj).bixialpha(ich, :).disp_stat();

    %% SAVE FOR EACH CHANNEL
    % save memory
    iEEG.subjinfo(isubj).hos(ich, :).compact_result();
    hos = iEEG.subjinfo(isubj).hos(ich, :);
    hos.coef = [];
    hos.sres = [];
    hos.bsres = [];

    bixialpha = iEEG.subjinfo(isubj).bixialpha(ich, :);

    savepath = fullfile(repo_root, 'result', 'iEEG', 'HoXiAlpha/', ...
        ['W', filesep, taskname, '_cpct', filesep], ...
        ['iEEG_', num2str(isubj), '_', num2str(ich), '.mat']);
    test_folder(savepath);
    save(savepath, "hos", "bixialpha", '-v7.3');

    % disp info
    iEEG.subjinfo(isubj).logging.t = toc;
    disp(['finish case ', num2str(isubj), '-', num2str(ich), ' ', iEEG.subjinfo(isubj).hos(ich, :).name, ', with ', num2str(iEEG.subjinfo(isubj).bixialpha(ich, :).logging.t), ' second'])
    disp(repmat(newline, 1, 3));

end

function [iEEG] = plot_save_result(iEEG, idx_subj_ch, repo_root, taskname, ids)

    for id = ids
        isubj = idx_subj_ch(id, 1);
        ich = idx_subj_ch(id, 2);
        try
            figname = [string(['Stage_', iEEG.subjinfo(isubj).hos(ich, :).info.stg, '_Region_', num2str(iEEG.subjinfo(isubj).hos(ich, :).info.rgid), ...
                '_Patient_', num2str(iEEG.subjinfo(isubj).hos(ich, :).info.ptid), '_Gender_', iEEG.subjinfo(isubj).hos(ich, :).info.gender, '_Age_', num2str(iEEG.subjinfo(isubj).hos(ich, :).info.age)]), ...
                string(['_CName_', iEEG.subjinfo(isubj).hos(ich, :).info.cname, '_CType_', iEEG.subjinfo(isubj).hos(ich, :).info.ctype, '_HemiS_', iEEG.subjinfo(isubj).hos(ich, :).info.hemis]), ...
                string(['_method_', iEEG.subjinfo(isubj).hos(ich, :).method, '_normlzt_', iEEG.subjinfo(isubj).hos(ich, :).normalization, '_Nw_', num2str(iEEG.subjinfo(isubj).hos(ich, :).NW), '_wind_', num2str(iEEG.subjinfo(isubj).hos(ich, :).window)]), ...
                string(['_', num2str(isubj), '_', num2str(ich)])];

            %% save result and then figure
            figpath = fullfile(repo_root, 'figure', 'iEEG', 'HoXiAlpha/', ...
                (iEEG.subjinfo(isubj).hos(ich, :).info.stg), filesep, ...
                taskname, filesep, 'bixialpha_', iEEG.subjinfo(isubj).ChName{ich}, '_xialpha_', num2str(isubj), '_', num2str(ich), '.jpg');

            if ~exist(figpath, 'file')
                f1 = figure(1);
                show(iEEG.subjinfo(isubj).bixialpha(ich, :).xialpha)
                figtitle = replace(strcat(figname, '_xialpha'), '_', '-');
                sgtitle(figtitle)

                test_folder(figpath);
                saveas(f1, figpath)
            end

            figpath = fullfile(repo_root, 'figure', 'iEEG', 'HoXiAlpha/', ...
                (iEEG.subjinfo(isubj).hos(ich, :).info.stg), filesep, ...
                taskname, filesep, 'bixialpha_', iEEG.subjinfo(isubj).ChName{ich}, '_stat_', num2str(isubj), '_', num2str(ich), '.jpg');

        catch ME
            disp(getReport(ME, 'extended', 'hyperlinks', 'on'))
            warning(['error in case ', iEEG.subjinfo(isubj).hos(ich, :).name])
        end

    end

end
