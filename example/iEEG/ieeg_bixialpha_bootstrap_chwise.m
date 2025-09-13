function ieeg_bixialpha_bootstrap_chwise(itask, npertask, taskname, iEEG)
currentTime = datetime('now');
disp(currentTime);
close all; %clc;

if ~exist('iEEG', 'var')
    filepath = './data/iEEG/iEEG1772FemTposHmInvScMleLog.mat';
    iEEG = load(filepath);
    iEEG = remove_useless_ieeg_fields(iEEG);
end

Nch_subj = cat_struct_fields(iEEG.subjinfo, 'NumChannel', false, true);
idx_subj_ch = gen_idx_subj_ch(Nch_subj);

Nid = length(idx_subj_ch);
[id, task] = task_take(taskname, Nid, itask, npertask);
disp(['[taskmng]: ', num2str(task.itask), 'th/', num2str(task.Ntask), ' task with ', num2str(task.npertask), ' ids: ', num2str(id)]);

startmatlabpool(feature('numcores'));
[iEEG] = batch_iEEG_hos_bixialpha(iEEG, idx_subj_ch, taskname, id);
disp(['[taskmng]:finish task']);

% plot result and save
if desktop('-inuse')
    plot_save_result(iEEG, idx_subj_ch, taskname, id);
end

currentTime = datetime('now');
disp(currentTime);
end

%%
function [iEEG] = batch_iEEG_hos_bixialpha(iEEG, idx_subj_ch, taskname, ids)
Fs = 200;

for id = ids
    isubj = idx_subj_ch(id, 1);
    ich = idx_subj_ch(id, 2);
    disp(['[taskmng]:start task ', num2str(isubj), '_', num2str(ich)]);
    get_freemem(true);
    savepath = ['./result/iEEG/HoXiAlpha/', ...
                    'W', filesep, taskname, '_cpct', filesep, ...
                    'iEEG_', num2str(isubj), '_', num2str(ich), '.mat'];

    if exist(savepath, "file")

        try
            % load(savepath, 'hos', 'nlng', 'bixialpha');
            load(savepath, 'hos', 'bixialpha');
            hos.compact_result();
            bixialpha.compact_result();
            iEEG.subjinfo(isubj).hos(ich, :) = hos;
            iEEG.subjinfo(isubj).bixialpha(ich, :) = bixialpha;
            % save(savepath, "hos", "nlng", "bixialpha")

        catch
            warning(['file crash ', num2str(isubj), '_', num2str(ich), ', try to recalculate'])
            iEEG = calc_iEEG_hos_bixialpha(iEEG, taskname, Fs, isubj, ich);
        end

    else
        warning(['lost file ', num2str(isubj), '_', num2str(ich), ', try to recalculate'])
        iEEG = calc_iEEG_hos_bixialpha(iEEG, taskname, Fs, isubj, ich);
    end

end

end

%%
function iEEG = calc_iEEG_hos_bixialpha(iEEG, taskname, Fs, isubj, ich)

try
    tic
    iEEG.subjinfo(isubj).hos(ich, :) = HOS(double(iEEG.subjinfo(isubj).Data_W(1:12000, ich)), Fs);
    iEEG.subjinfo(isubj).hos(ich, :).window = 300;
    iEEG.subjinfo(isubj).hos(ich, :).nfft = 300;
    iEEG.subjinfo(isubj).hos(ich, :).frange = [1.5, 45];

    iEEG.subjinfo(isubj).hos(ich, :).overlap = 0.75;
    iEEG.subjinfo(isubj).hos(ich, :).normalization = 'haubrich'; % haubrich skewness
    iEEG.subjinfo(isubj).hos(ich, :).NW = 1.5;
    iEEG.subjinfo(isubj).hos(ich, :).method = 'pmtm'; % fft pmtm
    iEEG.subjinfo(isubj).hos(ich, :).tf_return_res = false; %true false
    % iEEG.subjinfo(isubj).hos(ich, :).tf_return_seg = true; %true false
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
    iEEG.subjinfo(isubj).hos(ich, :).info.rgid = iEEG.subjinfo(isubj).ChannelRegionOrig(ich);
    iEEG.subjinfo(isubj).hos(ich, :).info.rgname = iEEG.RegionName(iEEG.subjinfo(isubj).hos(ich, :).info.rgid);
    iEEG.subjinfo(isubj).hos(ich, :).info.cpos = iEEG.subjinfo(isubj).ChannelPositionMid(ich, :);
    iEEG.subjinfo(isubj).hos(ich, :).info.gender = iEEG.subjinfo(isubj).Gender{1};
    iEEG.subjinfo(isubj).hos(ich, :).info.age = iEEG.subjinfo(isubj).AgeAtTimeOfStudy;

    ksall = 6;
    khall = 10;

    if isfield(iEEG.subjinfo(isubj), 'bixialpha') && isnumeric(iEEG.subjinfo(isubj).bixialpha)
        iEEG.subjinfo(isubj).bixialpha = BiXiAlphaBootstrap();
    end

    iEEG.subjinfo(isubj).bixialpha(ich, :) = BiXiAlphaBootstrap(iEEG.subjinfo(isubj).hos(ich, :), ksall, khall);
    iEEG.subjinfo(isubj).bixialpha(ich, :).multiple_seg = iEEG.subjinfo(isubj).hos(ich, :).tf_return_seg;
    iEEG.subjinfo(isubj).bixialpha(ich, :).verbose = false;
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

    iEEG.subjinfo(isubj).bixialpha(ich, :).compact_result();
    bixialpha = iEEG.subjinfo(isubj).bixialpha(ich, :);

    savepath = ['./result/iEEG/HoXiAlpha/', ...
                    iEEG.subjinfo(isubj).hos(ich).info.stg, filesep, taskname, '_cpct', filesep, ...
                    'iEEG_', num2str(isubj), '_', num2str(ich), '.mat'];
    test_folder(savepath);
    save(savepath, "hos", "bixialpha", '-v7.3');
    % save(savepath, "hos", "nlng", "bixialpha", '-v7.3');
    % disp info
    iEEG.subjinfo(isubj).logging.t = toc;
    disp(['finish case ', num2str(isubj), '-', num2str(ich), ' ', iEEG.subjinfo(isubj).hos(ich, :).name, ', with ', num2str(iEEG.subjinfo(isubj).bixialpha(ich, :).logging.t), ' second'])
    disp(repmat(newline, 1, 3));

catch ME
    disp(getReport(ME, 'extended', 'hyperlinks', 'on'))
    iEEG.subjinfo(isubj).bixialpha(ich, :).logging.log = ME;
    warning(['error in case ', iEEG.subjinfo(isubj).hos(ich, :).name])
    % print error in ME
    % warning(ME.message)
end

end

function [iEEG] = plot_save_result(iEEG, idx_subj_ch, taskname, ids)

for id = ids
    isubj = idx_subj_ch(id, 1);
    ich = idx_subj_ch(id, 2);
    % Nsubj=iEEG.Nsubj;
    try
        figname = [string(['Stage_', iEEG.subjinfo(isubj).hos(ich, :).info.stg, '_Region_', num2str(iEEG.subjinfo(isubj).hos(ich, :).info.rgid), ...
                               '_Patient_', num2str(iEEG.subjinfo(isubj).hos(ich, :).info.ptid), '_Gender_', iEEG.subjinfo(isubj).hos(ich, :).info.gender, '_Age_', num2str(iEEG.subjinfo(isubj).hos(ich, :).info.age)]), ...
            string(['_CName_', iEEG.subjinfo(isubj).hos(ich, :).info.cname, '_CType_', iEEG.subjinfo(isubj).hos(ich, :).info.ctype, '_HemiS_', iEEG.subjinfo(isubj).hos(ich, :).info.hemis]), ...
            string(['_method_', iEEG.subjinfo(isubj).hos(ich, :).method, '_normlzt_', iEEG.subjinfo(isubj).hos(ich, :).normalization, '_Nw_', num2str(iEEG.subjinfo(isubj).hos(ich, :).NW), '_wind_', num2str(iEEG.subjinfo(isubj).hos(ich, :).window)]), ...
            string(['_', num2str(isubj), '_', num2str(ich)])];

        %% save result and then figure
        figpath = ['./figure/iEEG/HoXiAlpha/', ...
                       (iEEG.subjinfo(isubj).hos(ich, :).info.stg), filesep, ...
                       taskname, filesep, 'bixialpha_', iEEG.subjinfo(isubj).ChName{ich}, '_xialpha_', num2str(isubj), '_', num2str(ich), '.jpg'];

        if ~exist(figpath, 'file')
            f1 = figure(1);
            show(iEEG.subjinfo(isubj).bixialpha(ich, :).xialpha)
            figtitle = replace(strcat(figname, '_xialpha'), '_', '-');
            sgtitle(figtitle)

            test_folder(figpath);
            saveas(f1, figpath)
        end

    catch ME
        disp(getReport(ME, 'extended', 'hyperlinks', 'on'))
        warning(['error in case ', iEEG.subjinfo(isubj).hos(ich, :).name])
    end

end

end
