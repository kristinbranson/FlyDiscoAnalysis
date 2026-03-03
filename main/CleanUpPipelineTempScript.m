experiment_names = {
'GMR_22E07_AE_01_TrpA_Rig1Plate10BowlC_20110608T130458'
'GMR_26D10_AE_01_TrpA_Rig2Plate17BowlA_20110819T112938'
'GMR_61D08_AD_01_TrpA_Rig1Plate15BowlA_20111006T160338'
'GMR_61D08_AD_01_TrpA_Rig1Plate15BowlB_20111006T160333'
'GMR_61G12_AD_01_TrpA_Rig2Plate17BowlD_20111006T160921'
'GMR_61H02_AD_01_TrpA_Rig1Plate15BowlA_20111006T162653'
'GMR_61H02_AD_01_TrpA_Rig1Plate15BowlB_20111006T162648'
'GMR_64B03_AE_01_TrpA_Rig1Plate15BowlA_20111007T133630'
'pBDPGAL4U_TrpA_Rig1Plate15BowlA_20111007T135905'
}'

experiment_names = cellfun(@(s) ['FlyBowl_',s],experiment_names,'UniformOutput',false)
data = SAGEListBowlExperiments('experiment_name',experiment_names,'checkflags',false)
[~,order] = sort({data.experiment_name})
data = data(order)
[{data.line_name}',{data.exp_datetime}',{data.automated_pf}',{data.automated_pf_category}']
setdiff(experiment_names,{data.experiment_name})