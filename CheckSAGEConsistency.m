%% check for consistency between experiment directory data and SAGE flattened view

%% set up path

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
end

%% parameters

% data sets to look at
datasets = {'score','data','histogram'};

% fields that should be available, per dataset
required_fields.score = {
  'experiment_id'
  'experiment_name'
  'experiment_protocol'
  'experimenter'
  'exp_datetime'
  'day_of_week'
  'automated_pf'
  'flag_aborted'
  'flag_redo'
  'flag_review'
  'manual_pf'
  'humidity'
  'notes_behavioral'
  'notes_technical'
  'temperature'
  'session_id'
  'session_name'
  'line_id'
  'line_name'
  'line_lab'
  'bowl'
  'camera'
  'computer'
  'harddrive'
  'apparatus_id'
  'cross_date'
  'effector'
  'environmental_chamber'
  'gender'
  'handler_sorting'
  'handler_starvation'
  'handling_protocol'
  'hours_sorted'
  'hours_starved'
  'num_flies'
  'num_flies_dead'
  'plate'
  'rearing_incubator'
  'rearing_protocol'
  'rig'
  'seconds_fliesloaded'
  'seconds_shiftflytemp'
  'top_plate'
  'num_flies_damaged'
  'handler_cross'
  'notes_keyword'
  'notes_curation'
  'cross_barcode'
  'flip_used'
  'flip_date'
  'room'
  'wish_list'
  'robot_stock_copy'
  'screen_type'
  'screen_reason'
  'data_capture_version'
  'analysis_protocol'
  'QuickStats_BackSubStats_meanNConnComps'
  'QuickStats_BackSubStats_minNConnComps'
  'QuickStats_BackSubStats_maxNConnComps'
  'QuickStats_BackSubStats_stdNConnComps'
  'QuickStats_BackSubStats_meanBlobArea'
  'QuickStats_BackSubStats_minBlobArea'
  'QuickStats_BackSubStats_maxBlobArea'
  'QuickStats_BackSubStats_stdBlobArea'
  'ufmf_diagnostics_summary_fracFramesWithFracFgPx0050000'
  'ufmf_diagnostics_summary_fracFramesWithFracFgPx0100000'
  'ufmf_diagnostics_summary_fracFramesWithFracFgPx0250000'
  'ufmf_diagnostics_summary_maxBandWidth'
  'ufmf_diagnostics_summary_maxComputeBackgroundTime'
  'ufmf_diagnostics_summary_maxComputeFrameTime'
  'ufmf_diagnostics_summary_maxComputeStatisticsTime'
  'ufmf_diagnostics_summary_maxFPS'
  'ufmf_diagnostics_summary_maxFrameSize'
  'ufmf_diagnostics_summary_maxMaxFilterError'
  'ufmf_diagnostics_summary_maxMaxPixelError'
  'ufmf_diagnostics_summary_maxMeanPixelError'
  'ufmf_diagnostics_summary_maxNBoxes'
  'ufmf_diagnostics_summary_maxNForegroundPx'
  'ufmf_diagnostics_summary_maxNPxWritten'
  'ufmf_diagnostics_summary_maxReweightCountsTime'
  'ufmf_diagnostics_summary_maxUpdateBackgroundTime'
  'ufmf_diagnostics_summary_maxWaitForFrameCaptureTime'
  'ufmf_diagnostics_summary_maxWriteFooterTime'
  'ufmf_diagnostics_summary_maxWriteFrameTime'
  'ufmf_diagnostics_summary_maxWriteHeaderTime'
  'ufmf_diagnostics_summary_maxWriteKeyFrameTime'
  'ufmf_diagnostics_summary_meanBandWidth'
  'ufmf_diagnostics_summary_meanCompressionRate'
  'ufmf_diagnostics_summary_meanComputeBackgroundTime'
  'ufmf_diagnostics_summary_meanComputeFrameTime'
  'ufmf_diagnostics_summary_meanComputeStatisticsTime'
  'ufmf_diagnostics_summary_meanFPS'
  'ufmf_diagnostics_summary_meanFrameSize'
  'ufmf_diagnostics_summary_meanMaxFilterError'
  'ufmf_diagnostics_summary_meanMaxPixelError'
  'ufmf_diagnostics_summary_meanMeanPixelError'
  'ufmf_diagnostics_summary_meanNBoxes'
  'ufmf_diagnostics_summary_meanNForegroundPx'
  'ufmf_diagnostics_summary_meanNPxWritten'
  'ufmf_diagnostics_summary_meanReweightCountsTime'
  'ufmf_diagnostics_summary_meanUpdateBackgroundTime'
  'ufmf_diagnostics_summary_meanWaitForFrameCaptureTime'
  'ufmf_diagnostics_summary_meanWriteFooterTime'
  'ufmf_diagnostics_summary_meanWriteFrameTime'
  'ufmf_diagnostics_summary_meanWriteHeaderTime'
  'ufmf_diagnostics_summary_meanWriteKeyFrameTime'
  'ufmf_diagnostics_summary_minFPS'
  'ufmf_diagnostics_summary_nComputeBackgroundCalls'
  'ufmf_diagnostics_summary_nComputeFrameCalls'
  'ufmf_diagnostics_summary_nComputeStatisticsCalls'
  'ufmf_diagnostics_summary_nFrames'
  'ufmf_diagnostics_summary_nFramesDroppedTotal'
  'ufmf_diagnostics_summary_nFramesNoBackSub'
  'ufmf_diagnostics_summary_nFramesUncompressed'
  'ufmf_diagnostics_summary_nReweightCountsCalls'
  'ufmf_diagnostics_summary_nUpdateBackgroundCalls'
  'ufmf_diagnostics_summary_nWaitForFrameCaptureCalls'
  'ufmf_diagnostics_summary_nWriteFooterCalls'
  'ufmf_diagnostics_summary_nWriteFrameCalls'
  'ufmf_diagnostics_summary_nWriteHeaderCalls'
  'ufmf_diagnostics_summary_nWriteKeyFrameCalls'
  'ufmf_diagnostics_summary_stdFPS'
  'ufmf_diagnostics_summary_stdFrameSize'
  'ufmf_diagnostics_summary_stdMaxFilterError'
  'ufmf_diagnostics_summary_stdMaxPixelError'
  'ufmf_diagnostics_summary_stdMeanPixelError'
  'ufmf_diagnostics_summary_stdNBoxes'
  'ufmf_diagnostics_summary_stdNForegroundPx'
  'ufmf_diagnostics_summary_stdNPxWritten'
  'ctrax_diagnostics_sum_nsplit'
  'ctrax_diagnostics_ndeaths_notfixed'
  'ctrax_diagnostics_nsmall_notfixed'
  'ctrax_diagnostics_nframes_analyzed'
  'ctrax_diagnostics_nlost_fixed'
  'ctrax_diagnostics_nmerged_fixed'
  'ctrax_diagnostics_nbirths_nohindsight'
  'ctrax_diagnostics_max_nsplit'
  'ctrax_diagnostics_nsplits_fixed'
  'ctrax_diagnostics_nsmall_lowerthresh'
  'ctrax_diagnostics_nlarge_notfixed'
  'ctrax_diagnostics_nsmall_merged'
  'ctrax_diagnostics_nlarge_split'
  'ctrax_diagnostics_nsmall_deleted'
  'ctrax_diagnostics_nbirths_notfixed'
  'ctrax_diagnostics_nlarge_ignored'
  'ctrax_diagnostics_ndeaths_nohindsight'
  'ctrax_diagnostics_nspurious_fixed'
  'ctrax_diagnostics_nhindsight_fixed'
  'registrationdata_offX'
  'registrationdata_offY'
  'registrationdata_offTheta'
  'registrationdata_scale'
  'registrationdata_bowlMarkerTheta'
  'registrationdata_featureStrengths'
  'registrationdata_circleCenterX'
  'registrationdata_circleCenterY'
  'registrationdata_circleRadius'
  'sexclassifier_diagnostics_mean_normhmmscore'
  'sexclassifier_diagnostics_median_normhmmscore'
  'sexclassifier_diagnostics_std_normhmmscore'
  'sexclassifier_diagnostics_min_normhmmscore'
  'sexclassifier_diagnostics_max_normhmmscore'
  'sexclassifier_diagnostics_mean_nswaps'
  'sexclassifier_diagnostics_median_nswaps'
  'sexclassifier_diagnostics_std_nswaps'
  'sexclassifier_diagnostics_min_nswaps'
  'sexclassifier_diagnostics_max_nswaps'
  'sexclassifier_diagnostics_mean_meanabsdev'
  'sexclassifier_diagnostics_median_meanabsdev'
  'sexclassifier_diagnostics_std_meanabsdev'
  'sexclassifier_diagnostics_min_meanabsdev'
  'sexclassifier_diagnostics_max_meanabsdev'
  'sexclassifier_diagnostics_mean_nfemales'
  'sexclassifier_diagnostics_median_nfemales'
  'sexclassifier_diagnostics_std_nfemales'
  'sexclassifier_diagnostics_min_nfemales'
  'sexclassifier_diagnostics_max_nfemales'
  'sexclassifier_diagnostics_mean_nmales'
  'sexclassifier_diagnostics_median_nmales'
  'sexclassifier_diagnostics_std_nmales'
  'sexclassifier_diagnostics_min_nmales'
  'sexclassifier_diagnostics_max_nmales'
  'sexclassifier_diagnostics_mean_nflies'
  'sexclassifier_diagnostics_median_nflies'
  'sexclassifier_diagnostics_std_nflies'
  'sexclassifier_diagnostics_min_nflies'
  'sexclassifier_diagnostics_max_nflies'
  'temperature_diagnostics_mean'
  'temperature_diagnostics_max'
  'temperature_diagnostics_maxdiff'
  'temperature_diagnostics_nreadings'
  'temperature_diagnostics_std'
  'bkgd_diagnostics_histmode_bkgdcenter'
  'bkgd_diagnostics_mean_bkgdcenter'
  'bkgd_diagnostics_std_bkgdcenter'
  'bkgd_diagnostics_histmode_bkgddev'
  'bkgd_diagnostics_mean_bkgddev'
  'bkgd_diagnostics_std_bkgddev'
  'bkgd_diagnostics_histmode_alwaysbkgd'
  'bkgd_diagnostics_mean_alwaysbkgd'
  'bkgd_diagnostics_std_alwaysbkgd'
  'bkgd_diagnostics_histmode_bkgdcenter_llr'
  'bkgd_diagnostics_mean_bkgdcenter_llr'
  'bkgd_diagnostics_std_bkgdcenter_llr'
  'bkgd_diagnostics_histmode_imfore'
  'bkgd_diagnostics_mean_imfore'
  'bkgd_diagnostics_std_imfore'
  'bkgd_diagnostics_histmode_diffim'
  'bkgd_diagnostics_mean_diffim'
  'bkgd_diagnostics_std_diffim'
  'bkgd_diagnostics_histmode_llrfore'
  'bkgd_diagnostics_mean_llrfore'
  'bkgd_diagnostics_std_llrfore'
  'bkgd_diagnostics_histmode_minllrperfly'
  'bkgd_diagnostics_mean_minllrperfly'
  'bkgd_diagnostics_histmode_meanllrperfly'
  'bkgd_diagnostics_mean_meanllrperfly'
  'bkgd_diagnostics_histmode_maxllrperfly'
  'bkgd_diagnostics_mean_maxllrperfly'
  'bias_diagnostics_max_sumfracsmooth'
  'bias_diagnostics_min_sumfracsmooth'
  'bias_diagnostics_argmaxangle_sumfracsmooth'
  'bias_diagnostics_argminangle_sumfracsmooth'
  'bias_diagnostics_maxratio_sumfracsmooth'
  'bias_diagnostics_diffargextremaangle_sumfracsmooth'
  'bias_diagnostics_max_sumfracallsmooth'
  'bias_diagnostics_min_sumfracallsmooth'
  'bias_diagnostics_argmaxangle_sumfracallsmooth'
  'bias_diagnostics_argminangle_sumfracallsmooth'
  'bias_diagnostics_maxratio_sumfracallsmooth'
  'bias_diagnostics_diffargextremaangle_sumfracallsmooth'
  'bias_diagnostics_max_maxfracsmooth'
  'bias_diagnostics_min_maxfracsmooth'
  'bias_diagnostics_argmaxangle_maxfracsmooth'
  'bias_diagnostics_argminangle_minfracsmooth'
  'bias_diagnostics_maxdiff_maxfracsmooth'
  'bias_diagnostics_diffargextremaangle_maxfracsmooth'
  'bias_diagnostics_max_maxfracallsmooth'
  'bias_diagnostics_min_maxfracallsmooth'
  'bias_diagnostics_argmaxangle_maxfracallsmooth'
  'bias_diagnostics_argminangle_minfracallsmooth'
  'bias_diagnostics_maxdiff_maxfracallsmooth'
  'bias_diagnostics_diffargextremaangle_maxfracallsmooth'
  'registrationdata_end_frame'
  'registrationdata_seconds_crop_end'
  'registrationdata_seconds_crop_start'
  'registrationdata_start_frame'
  'sexclassifier_diagnostics_classifier_mu_area_female'
  'sexclassifier_diagnostics_classifier_mu_area_male'
  'sexclassifier_diagnostics_classifier_var_area_female'
  'sexclassifier_diagnostics_classifier_var_area_male'
  'sexclassifier_diagnostics_classifier_loglik'
  'sexclassifier_diagnostics_classifier_niters'
  'rig_bowl'
  'line__effector'
  };

required_fields.data = {
  'experiment_id'
  'experiment_name'
  'experiment_protocol'
  'experimenter'
  'exp_datetime'
  'file_system_path'
  'day_of_week'
  'automated_pf'
  'flag_aborted'
  'flag_redo'
  'flag_review'
  'manual_pf'
  'humidity'
  'notes_behavioral'
  'notes_technical'
  'temperature'
  'session_id'
  'session_name'
  'line_id'
  'line_name'
  'line_lab'
  'bowl'
  'camera'
  'computer'
  'harddrive'
  'apparatus_id'
  'cross_date'
  'effector'
  'environmental_chamber'
  'gender'
  'handler_sorting'
  'handler_starvation'
  'handling_protocol'
  'hours_sorted'
  'hours_starved'
  'num_flies'
  'num_flies_dead'
  'plate'
  'rearing_incubator'
  'rearing_protocol'
  'rig'
  'seconds_fliesloaded'
  'seconds_shiftflytemp'
  'top_plate'
  'num_flies_damaged'
  'handler_cross'
  'notes_keyword'
  'notes_curation'
  'cross_barcode'
  'flip_used'
  'flip_date'
  'room'
  'wish_list'
  'robot_stock_copy'
  'screen_type'
  'screen_reason'
  'data_capture_version'
  'analysis_protocol'
  'QuickStats_BackSubStats_meanNConnComps'
  'QuickStats_BackSubStats_minNConnComps'
  'QuickStats_BackSubStats_maxNConnComps'
  'QuickStats_BackSubStats_stdNConnComps'
  'QuickStats_BackSubStats_meanBlobArea'
  'QuickStats_BackSubStats_minBlobArea'
  'QuickStats_BackSubStats_maxBlobArea'
  'QuickStats_BackSubStats_stdBlobArea'
  'ufmf_diagnostics_summary_fracFramesWithFracFgPx0050000'
  'ufmf_diagnostics_summary_fracFramesWithFracFgPx0100000'
  'ufmf_diagnostics_summary_fracFramesWithFracFgPx0250000'
  'ufmf_diagnostics_summary_maxBandWidth'
  'ufmf_diagnostics_summary_maxComputeBackgroundTime'
  'ufmf_diagnostics_summary_maxComputeFrameTime'
  'ufmf_diagnostics_summary_maxComputeStatisticsTime'
  'ufmf_diagnostics_summary_maxFPS'
  'ufmf_diagnostics_summary_maxFrameSize'
  'ufmf_diagnostics_summary_maxMaxFilterError'
  'ufmf_diagnostics_summary_maxMaxPixelError'
  'ufmf_diagnostics_summary_maxMeanPixelError'
  'ufmf_diagnostics_summary_maxNBoxes'
  'ufmf_diagnostics_summary_maxNForegroundPx'
  'ufmf_diagnostics_summary_maxNPxWritten'
  'ufmf_diagnostics_summary_maxReweightCountsTime'
  'ufmf_diagnostics_summary_maxUpdateBackgroundTime'
  'ufmf_diagnostics_summary_maxWaitForFrameCaptureTime'
  'ufmf_diagnostics_summary_maxWriteFooterTime'
  'ufmf_diagnostics_summary_maxWriteFrameTime'
  'ufmf_diagnostics_summary_maxWriteHeaderTime'
  'ufmf_diagnostics_summary_maxWriteKeyFrameTime'
  'ufmf_diagnostics_summary_meanBandWidth'
  'ufmf_diagnostics_summary_meanCompressionRate'
  'ufmf_diagnostics_summary_meanComputeBackgroundTime'
  'ufmf_diagnostics_summary_meanComputeFrameTime'
  'ufmf_diagnostics_summary_meanComputeStatisticsTime'
  'ufmf_diagnostics_summary_meanFPS'
  'ufmf_diagnostics_summary_meanFrameSize'
  'ufmf_diagnostics_summary_meanMaxFilterError'
  'ufmf_diagnostics_summary_meanMaxPixelError'
  'ufmf_diagnostics_summary_meanMeanPixelError'
  'ufmf_diagnostics_summary_meanNBoxes'
  'ufmf_diagnostics_summary_meanNForegroundPx'
  'ufmf_diagnostics_summary_meanNPxWritten'
  'ufmf_diagnostics_summary_meanReweightCountsTime'
  'ufmf_diagnostics_summary_meanUpdateBackgroundTime'
  'ufmf_diagnostics_summary_meanWaitForFrameCaptureTime'
  'ufmf_diagnostics_summary_meanWriteFooterTime'
  'ufmf_diagnostics_summary_meanWriteFrameTime'
  'ufmf_diagnostics_summary_meanWriteHeaderTime'
  'ufmf_diagnostics_summary_meanWriteKeyFrameTime'
  'ufmf_diagnostics_summary_minFPS'
  'ufmf_diagnostics_summary_nComputeBackgroundCalls'
  'ufmf_diagnostics_summary_nComputeFrameCalls'
  'ufmf_diagnostics_summary_nComputeStatisticsCalls'
  'ufmf_diagnostics_summary_nFrames'
  'ufmf_diagnostics_summary_nFramesDroppedTotal'
  'ufmf_diagnostics_summary_nFramesNoBackSub'
  'ufmf_diagnostics_summary_nFramesUncompressed'
  'ufmf_diagnostics_summary_nReweightCountsCalls'
  'ufmf_diagnostics_summary_nUpdateBackgroundCalls'
  'ufmf_diagnostics_summary_nWaitForFrameCaptureCalls'
  'ufmf_diagnostics_summary_nWriteFooterCalls'
  'ufmf_diagnostics_summary_nWriteFrameCalls'
  'ufmf_diagnostics_summary_nWriteHeaderCalls'
  'ufmf_diagnostics_summary_nWriteKeyFrameCalls'
  'ufmf_diagnostics_summary_stdFPS'
  'ufmf_diagnostics_summary_stdFrameSize'
  'ufmf_diagnostics_summary_stdMaxFilterError'
  'ufmf_diagnostics_summary_stdMaxPixelError'
  'ufmf_diagnostics_summary_stdMeanPixelError'
  'ufmf_diagnostics_summary_stdNBoxes'
  'ufmf_diagnostics_summary_stdNForegroundPx'
  'ufmf_diagnostics_summary_stdNPxWritten'
  'ctrax_diagnostics_sum_nsplit'
  'ctrax_diagnostics_ndeaths_notfixed'
  'ctrax_diagnostics_nsmall_notfixed'
  'ctrax_diagnostics_nframes_analyzed'
  'ctrax_diagnostics_nlost_fixed'
  'ctrax_diagnostics_nmerged_fixed'
  'ctrax_diagnostics_nbirths_nohindsight'
  'ctrax_diagnostics_max_nsplit'
  'ctrax_diagnostics_nsplits_fixed'
  'ctrax_diagnostics_nsmall_lowerthresh'
  'ctrax_diagnostics_nlarge_notfixed'
  'ctrax_diagnostics_nsmall_merged'
  'ctrax_diagnostics_nlarge_split'
  'ctrax_diagnostics_nsmall_deleted'
  'ctrax_diagnostics_nbirths_notfixed'
  'ctrax_diagnostics_nlarge_ignored'
  'ctrax_diagnostics_ndeaths_nohindsight'
  'ctrax_diagnostics_nspurious_fixed'
  'ctrax_diagnostics_nhindsight_fixed'
  'registrationdata_offX'
  'registrationdata_offY'
  'registrationdata_offTheta'
  'registrationdata_scale'
  'registrationdata_bowlMarkerTheta'
  'registrationdata_featureStrengths'
  'registrationdata_circleCenterX'
  'registrationdata_circleCenterY'
  'registrationdata_circleRadius'
  'sexclassifier_diagnostics_mean_normhmmscore'
  'sexclassifier_diagnostics_median_normhmmscore'
  'sexclassifier_diagnostics_std_normhmmscore'
  'sexclassifier_diagnostics_min_normhmmscore'
  'sexclassifier_diagnostics_max_normhmmscore'
  'sexclassifier_diagnostics_mean_nswaps'
  'sexclassifier_diagnostics_median_nswaps'
  'sexclassifier_diagnostics_std_nswaps'
  'sexclassifier_diagnostics_min_nswaps'
  'sexclassifier_diagnostics_max_nswaps'
  'sexclassifier_diagnostics_mean_meanabsdev'
  'sexclassifier_diagnostics_median_meanabsdev'
  'sexclassifier_diagnostics_std_meanabsdev'
  'sexclassifier_diagnostics_min_meanabsdev'
  'sexclassifier_diagnostics_max_meanabsdev'
  'sexclassifier_diagnostics_mean_nfemales'
  'sexclassifier_diagnostics_median_nfemales'
  'sexclassifier_diagnostics_std_nfemales'
  'sexclassifier_diagnostics_min_nfemales'
  'sexclassifier_diagnostics_max_nfemales'
  'sexclassifier_diagnostics_mean_nmales'
  'sexclassifier_diagnostics_median_nmales'
  'sexclassifier_diagnostics_std_nmales'
  'sexclassifier_diagnostics_min_nmales'
  'sexclassifier_diagnostics_max_nmales'
  'sexclassifier_diagnostics_mean_nflies'
  'sexclassifier_diagnostics_median_nflies'
  'sexclassifier_diagnostics_std_nflies'
  'sexclassifier_diagnostics_min_nflies'
  'sexclassifier_diagnostics_max_nflies'
  'temperature_diagnostics_mean'
  'temperature_diagnostics_max'
  'temperature_diagnostics_maxdiff'
  'temperature_diagnostics_nreadings'
  'temperature_diagnostics_std'
  'bkgd_diagnostics_histmode_bkgdcenter'
  'bkgd_diagnostics_mean_bkgdcenter'
  'bkgd_diagnostics_std_bkgdcenter'
  'bkgd_diagnostics_histmode_bkgddev'
  'bkgd_diagnostics_mean_bkgddev'
  'bkgd_diagnostics_std_bkgddev'
  'bkgd_diagnostics_histmode_alwaysbkgd'
  'bkgd_diagnostics_mean_alwaysbkgd'
  'bkgd_diagnostics_std_alwaysbkgd'
  'bkgd_diagnostics_histmode_bkgdcenter_llr'
  'bkgd_diagnostics_mean_bkgdcenter_llr'
  'bkgd_diagnostics_std_bkgdcenter_llr'
  'bkgd_diagnostics_histmode_imfore'
  'bkgd_diagnostics_mean_imfore'
  'bkgd_diagnostics_std_imfore'
  'bkgd_diagnostics_histmode_diffim'
  'bkgd_diagnostics_mean_diffim'
  'bkgd_diagnostics_std_diffim'
  'bkgd_diagnostics_histmode_llrfore'
  'bkgd_diagnostics_mean_llrfore'
  'bkgd_diagnostics_std_llrfore'
  'bkgd_diagnostics_histmode_minllrperfly'
  'bkgd_diagnostics_mean_minllrperfly'
  'bkgd_diagnostics_histmode_meanllrperfly'
  'bkgd_diagnostics_mean_meanllrperfly'
  'bkgd_diagnostics_histmode_maxllrperfly'
  'bkgd_diagnostics_mean_maxllrperfly'
  'bias_diagnostics_max_sumfracsmooth'
  'bias_diagnostics_min_sumfracsmooth'
  'bias_diagnostics_argmaxangle_sumfracsmooth'
  'bias_diagnostics_argminangle_sumfracsmooth'
  'bias_diagnostics_maxratio_sumfracsmooth'
  'bias_diagnostics_diffargextremaangle_sumfracsmooth'
  'bias_diagnostics_max_sumfracallsmooth'
  'bias_diagnostics_min_sumfracallsmooth'
  'bias_diagnostics_argmaxangle_sumfracallsmooth'
  'bias_diagnostics_argminangle_sumfracallsmooth'
  'bias_diagnostics_maxratio_sumfracallsmooth'
  'bias_diagnostics_diffargextremaangle_sumfracallsmooth'
  'bias_diagnostics_max_maxfracsmooth'
  'bias_diagnostics_min_maxfracsmooth'
  'bias_diagnostics_argmaxangle_maxfracsmooth'
  'bias_diagnostics_argminangle_minfracsmooth'
  'bias_diagnostics_maxdiff_maxfracsmooth'
  'bias_diagnostics_diffargextremaangle_maxfracsmooth'
  'bias_diagnostics_max_maxfracallsmooth'
  'bias_diagnostics_min_maxfracallsmooth'
  'bias_diagnostics_argmaxangle_maxfracallsmooth'
  'bias_diagnostics_argminangle_minfracallsmooth'
  'bias_diagnostics_maxdiff_maxfracallsmooth'
  'bias_diagnostics_diffargextremaangle_maxfracallsmooth'
  'registrationdata_end_frame'
  'registrationdata_seconds_crop_end'
  'registrationdata_seconds_crop_start'
  'registrationdata_start_frame'
  'sexclassifier_diagnostics_classifier_mu_area_female'
  'sexclassifier_diagnostics_classifier_mu_area_male'
  'sexclassifier_diagnostics_classifier_var_area_female'
  'sexclassifier_diagnostics_classifier_var_area_male'
  'sexclassifier_diagnostics_classifier_loglik'
  'sexclassifier_diagnostics_classifier_niters'
  'QuickStats_BkgdIntensityHist_ctrs'
  'QuickStats_BkgdIntensityHist_frac'
  'QuickStats_BkgdScanLine_theta'
  'QuickStats_BkgdScanLine_intensities_1'
  'QuickStats_BkgdScanLine_intensities_2'
  'QuickStats_BkgdScanLine_intensities_3'
  'QuickStats_BkgdScanLine_intensities_4'
  'ufmf_diagnostics_stream_FPS'
  'ufmf_diagnostics_stream_bytes'
  'ufmf_diagnostics_stream_computeBackgroundTime'
  'ufmf_diagnostics_stream_computeFrameTime'
  'ufmf_diagnostics_stream_computeStatisticsTime'
  'ufmf_diagnostics_stream_frame'
  'ufmf_diagnostics_stream_isCompressed'
  'ufmf_diagnostics_stream_maxFilterError'
  'ufmf_diagnostics_stream_maxPixelError'
  'ufmf_diagnostics_stream_meanPixelError'
  'ufmf_diagnostics_stream_nBoxes'
  'ufmf_diagnostics_stream_nForegroundPx'
  'ufmf_diagnostics_stream_nFramesBuffered'
  'ufmf_diagnostics_stream_nFramesDropped'
  'ufmf_diagnostics_stream_nPxWritten'
  'ufmf_diagnostics_stream_reweightCountsTime'
  'ufmf_diagnostics_stream_timestamp'
  'ufmf_diagnostics_stream_updateBackgroundTime'
  'ufmf_diagnostics_stream_waitForFrameCaptureTime'
  'ufmf_diagnostics_stream_writeFooterTime'
  'ufmf_diagnostics_stream_writeFrameTime'
  'ufmf_diagnostics_stream_writeHeaderTime'
  'ufmf_diagnostics_stream_writeKeyFrameTime'
  'temperature_stream'
  'bkgd_diagnostics_prctiles_bkgdcenter'
  'bkgd_diagnostics_frac_bkgdcenter'
  'bkgd_diagnostics_prctiles_bkgddev'
  'bkgd_diagnostics_frac_bkgddev'
  'bkgd_diagnostics_prctiles_alwaysbkgd'
  'bkgd_diagnostics_frac_alwaysbkgd'
  'bkgd_diagnostics_prctiles_bkgdcenter_llr'
  'bkgd_diagnostics_frac_bkgdcenter_llr'
  'bkgd_diagnostics_prctiles_imfore'
  'bkgd_diagnostics_frac_imfore'
  'bkgd_diagnostics_prctiles_diffim'
  'bkgd_diagnostics_frac_diffim'
  'bkgd_diagnostics_prctiles_llrfore'
  'bkgd_diagnostics_frac_llrfore'
  'bkgd_diagnostics_frac_minllrperfly'
  'bkgd_diagnostics_frac_meanllrperfly'
  'bkgd_diagnostics_frac_maxllrperfly'
  'bias_diagnostics_frac'
  'bias_diagnostics_fracsmooth'
  'bias_diagnostics_fracall'
  'bias_diagnostics_fracallsmooth'
  'bias_diagnostics_sumfracsmooth'
  'bias_diagnostics_sumfracallsmooth'
  'bias_diagnostics_maxfracsmooth'
  'bias_diagnostics_maxfracallsmooth'
  'bias_diagnostics_argmaxr_fracsmooth'
  'bias_diagnostics_argmaxr_fracallsmooth'
  'rig_bowl'
  'line__effector'
  'stats_perframe_a_mm'
  'stats_perframe_absangle2wall'
  'stats_perframe_absanglefrom1to2_anglesub'
  'stats_perframe_absanglefrom1to2_nose2ell'
  'stats_perframe_absdangle2wall'
  'stats_perframe_absdtheta'
  'stats_perframe_absdv_cor'
  'stats_perframe_absphidiff_anglesub'
  'stats_perframe_absphidiff_nose2ell'
  'stats_perframe_absthetadiff_anglesub'
  'stats_perframe_absthetadiff_nose2ell'
  'stats_perframe_absyaw'
  'stats_perframe_anglesub'
  'stats_perframe_areasmooth'
  'stats_perframe_b_mm'
  'stats_perframe_corfrac_maj'
  'stats_perframe_dangle2wall'
  'stats_perframe_danglesub'
  'stats_perframe_dcenter'
  'stats_perframe_ddcenter'
  'stats_perframe_ddell2nose'
  'stats_perframe_ddist2wall'
  'stats_perframe_ddnose2ell'
  'stats_perframe_dell2nose'
  'stats_perframe_dist2wall'
  'stats_perframe_dnose2ell'
  'stats_perframe_dtheta'
  'stats_perframe_du_cor'
  'stats_perframe_du_ctr'
  'stats_perframe_dv_cor'
  'stats_perframe_dv_ctr'
  'stats_perframe_flipdv_cor'
  'stats_perframe_magveldiff_anglesub'
  'stats_perframe_magveldiff_nose2ell'
  'stats_perframe_phisideways'
  'stats_perframe_velmag'
  'stats_perframe_velmag_ctr'
  'stats_perframe_veltoward_anglesub'
  'stats_perframe_veltoward_nose2ell'
  'stats_perframe_x_mm'
  'stats_perframe_y_mm'
  'stats_perframe_yaw'
  };

required_fields.histogram = {
      'experiment_id'
    'experiment_name'
    'experiment_protocol'
    'experimenter'
    'exp_datetime'
    'file_system_path'
    'day_of_week'
    'automated_pf'
    'flag_aborted'
    'flag_redo'
    'flag_review'
    'manual_pf'
    'humidity'
    'notes_behavioral'
    'notes_technical'
    'temperature'
    'session_id'
    'session_name'
    'line_id'
    'line_name'
    'line_lab'
    'bowl'
    'camera'
    'computer'
    'harddrive'
    'apparatus_id'
    'cross_date'
    'effector'
    'environmental_chamber'
    'gender'
    'handler_sorting'
    'handler_starvation'
    'handling_protocol'
    'hours_sorted'
    'hours_starved'
    'num_flies'
    'num_flies_dead'
    'plate'
    'rearing_incubator'
    'rearing_protocol'
    'rig'
    'seconds_fliesloaded'
    'seconds_shiftflytemp'
    'top_plate'
    'num_flies_damaged'
    'handler_cross'
    'notes_keyword'
    'notes_curation'
    'cross_barcode'
    'flip_used'
    'flip_date'
    'room'
    'wish_list'
    'robot_stock_copy'
    'screen_type'
    'screen_reason'
    'data_capture_version'
    'analysis_protocol'
    'rig_bowl'
    'line__effector'
    'hist_perframe_a_mm'
    'hist_perframe_absangle2wall'
    'hist_perframe_absanglefrom1to2_anglesub'
    'hist_perframe_absanglefrom1to2_nose2ell'
    'hist_perframe_absdangle2wall'
    'hist_perframe_absdtheta'
    'hist_perframe_absdv_cor'
    'hist_perframe_absphidiff_anglesub'
    'hist_perframe_absphidiff_nose2ell'
    'hist_perframe_absthetadiff_anglesub'
    'hist_perframe_absthetadiff_nose2ell'
    'hist_perframe_absyaw'
    'hist_perframe_anglesub'
    'hist_perframe_areasmooth'
    'hist_perframe_arena_angle'
    'hist_perframe_b_mm'
    'hist_perframe_corfrac_maj'
    'hist_perframe_dangle2wall'
    'hist_perframe_danglesub'
    'hist_perframe_dcenter'
    'hist_perframe_ddcenter'
    'hist_perframe_ddell2nose'
    'hist_perframe_ddist2wall'
    'hist_perframe_ddnose2ell'
    'hist_perframe_dell2nose'
    'hist_perframe_dist2wall'
    'hist_perframe_dnose2ell'
    'hist_perframe_dtheta'
    'hist_perframe_du_cor'
    'hist_perframe_du_ctr'
    'hist_perframe_dv_cor'
    'hist_perframe_dv_ctr'
    'hist_perframe_flipdv_cor'
    'hist_perframe_magveldiff_anglesub'
    'hist_perframe_magveldiff_nose2ell'
    'hist_perframe_phi'
    'hist_perframe_phisideways'
    'hist_perframe_theta_mm'
    'hist_perframe_velmag'
    'hist_perframe_velmag_ctr'
    'hist_perframe_veltoward_anglesub'
    'hist_perframe_veltoward_nose2ell'
    'hist_perframe_x_mm'
    'hist_perframe_y_mm'
    'hist_perframe_yaw'
};

ntries = 3;

%% find all experiment directories that are successes

allexpdirs = dir(rootdatadir);
expdirs = {};
for expi = 1:numel(allexpdirs),

  expdir = fullfile(rootdatadir,allexpdirs(expi).name);
  if exist(fullfile(expdir,'ABORTED'),'file'),
    fprintf('Skipping experiment %s, ABORTED.\n',allexpdirs(expi).name);
    continue;
  elseif ~exist(fullfile(expdir,'SUCCESS'),'file'),
    fprintf('Skipping experiment %s, no SUCCESS file.\n',allexpdirs(expi).name);
    continue;
  end
  expdirs{end+1} = expdir; %#ok<SAGROW>
end
nexps = numel(expdirs);

%% query SAGE for each experiment

missingfns = struct;
extrafns = struct;
missingexps = struct;
failures = struct;
warningexps = struct;
for dataseti = 1:numel(datasets),
  missingfns.(datasets{dataseti}) = cell(1,nexps);
  extrafns.(datasets{dataseti}) = cell(1,nexps);
  missingexps.(datasets{dataseti}) = {};
  failures.(datasets{dataseti}) = {};
  warningexps.(datasets{dataseti}) = {};
end

%% pull each experiment

for expi = expi:nexps,

  if mod(expi,10) == 0,
    fprintf('Experiment %d / %d\n',expi,nexps);
  end
  
  [~,experiment_name] = fileparts(expdirs{expi});
  experiment_name = ['FlyBowl_',experiment_name]; %#ok<AGROW>
  
  for dataseti = 1:numel(datasets),
    dataset = datasets{dataseti};
    isfailure = false;
    for ntry = 1:ntries,
      try
        [tmp,~,iswarning] = SAGEGetBowlData('experiment_name',experiment_name,'dataset',dataset,'removemissingdata',false);
        break;
      catch ME,
        warning('try %d: SAGEGetBowlData error: %s',ntry,getReport(ME));
        if ntry == ntries,
          failures.(dataset){end+1} = experiment_name;
          isfailure = true;
        end
      end
    end
    if isfailure,
      continue;
    elseif iswarning,
      warningexps.(dataset){end+1} = experiment_name;
      continue;
    elseif isempty(tmp),
      missingexps.(dataset){end+1} = experiment_name;
      continue;
    end
    fns = fieldnames(tmp);
    missingfns.(dataset){expi} = setdiff(required_fields.(dataset),fns);
    extrafns.(dataset){expi} = setdiff(fns,required_fields.(dataset));
  end
  
  if mod(expi,100) == 0,
    save ConsistencyCheck.mat expi missingfns extrafns missingexps failures expdirs
  end
  
end

%% print results

filename = sprintf('SAGEConsistencyCheck%s.txt',datestr(now,'yyyymmdd'));
fid = fopen(filename,'w');

for dataseti = 1:numel(datasets),
  dataset = datasets{dataseti};
  for expi = find(~cellfun(@isempty,missingfns.(dataset))),
    [~,experiment_name] = fileparts(expdirs{expi});
    fprintf(fid,'%s %s missing fields:\n',experiment_name,dataset);
    fprintf(fid,'%s\n',missingfns.(dataset){expi}{:});
  end
end

for dataseti = 1:numel(datasets),
  
  dataset = datasets{dataseti};
  tmp = unique([extrafns.(dataset){~cellfun(@isempty,extrafns.(dataset))}]);
  if ~isempty(tmp),
    fprintf(fid,'Extra fields for %s:\n',dataset);
    fprintf(fid,'%s\n',tmp{:});
  end
  
end

for dataseti = 1:numel(datasets),
  dataset = datasets{dataseti};
  if ~isempty(missingexps.(dataset)),
    fprintf(fid,'Missing experiments for data set %s:\n',dataset);
    fprintf(fid,'%s\n',missingexps.(dataset){:});
  end
  if ~isempty(failures.(dataset)),
    fprintf(fid,'Errors for data set %s, experiments:\n',dataset);
    fprintf(fid,'%s\n',failures.(dataset){:});
  end
  if ~isempty(warningexps.(dataset)),
    fprintf(fid,'Warnings for data set %s, experiments:\n',dataset);
    fprintf(fid,'%s\n',warningexps.(dataset){:});
  end
end

fclose(fid);

%% fixing bugs above without rerunning everything
% 
% missingexps = struct;
% failures = struct;
% warningexps = struct;
% for dataseti = 1:numel(datasets),
%   missingexps.(datasets{dataseti}) = {};
%   failures.(datasets{dataseti}) = {};
%   warningexps.(datasets{dataseti}) = {};
% end
% 
% for dataseti = 1:numel(datasets),
%   dataset = datasets{dataseti};
%   for i = 1:numel(missingexps0.(dataset)),
%     fprintf('%s: %d/%d...\n',dataset,i,numel(missingexps0.(dataset)));
%     experiment_name = missingexps0.(dataset){i};
%     isfailure = false;
%     for ntry = 1:ntries,
%       try
%         [tmp,~,iswarning] = SAGEGetBowlData('experiment_name',experiment_name,'dataset',dataset,'removemissingdata',false);
%         break;
%       catch ME,
%         warning('try %d: SAGEGetBowlData error: %s',ntry,getReport(ME));
%         if ntry == ntries,
%           failures.(dataset){end+1} = experiment_name;
%           isfailure = true;
%           fprintf('Failed to load %s\n',experiment_name);
%         end
%       end
%     end
%     if isfailure,
%     elseif iswarning,
%       warningexps.(dataset){end+1} = experiment_name;
%       fprintf('Warning generated for experiment %s\n',experiment_name);
%     elseif isempty(tmp),
%       missingexps.(dataset){end+1} = experiment_name;
%       fprintf('No data for experiment %s\n',experiment_name);
%     else
%       fprintf('Experiment %s seems okay afterall\n',experiment_name);
%     end
%   end
% end