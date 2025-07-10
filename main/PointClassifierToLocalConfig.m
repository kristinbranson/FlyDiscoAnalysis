function PointClassifierToLocalConfig(classifierfile,configfile)

tmp = load(classifierfile);
tmp.configfilename = configfile;
save(classifierfile,'-struct','tmp');