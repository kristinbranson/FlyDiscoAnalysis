rootdatadir = '/groups/branson/bransonlab/flydisco_data';
outrootdatadir = '/groups/branson/home/robiea/Projects_data/JAABA/Data_FlyDisco';

expnamelist = {'VNC2_YNA_K_162984_RigB_20231205T114631'
    'VNC2_YNA_K_162984_RigA_20231206T125420'
    'VNC2_YNA_K_162984_RigA_20231205T114519'
    'VNC2_YNA_K_162984_RigB_20231206T125502'
    'VNC2_JRC_SS87693_RigB_20231205T112525'
    'VNC2_JRC_SS89300_RigA_20231205T120627'
    'VNC2_JRC_SS87693_RigD_20231206T110427'
    'VNC2_JRC_SS89300_RigC_20231206T130806'};

expnamelist = {'VNC2_JRC_SS72061_RigB_20230614T120606'
    'VNC2_JRC_SS74266_RigD_20230615T115903'}

for i = 1:numel(expnamelist)
indir = fullfile(rootdatadir,expnamelist{i});
outdir = fullfile(outrootdatadir,expnamelist{i});
copyfile(indir,outdir)

end
