% ScriptCheckExperiments20120726

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlDataCapture;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/starvation;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/age;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
rootdatadir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120328_non_olympiad_albins/results';
analysis_protocol = 'current_non_olympiad_albins';
npar = 8;
outexpdirfilename = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120328_non_olympiad_albins/expdirs_albins_20120911.txt';
expdirs = importdata(outexpdirfilename);
classifierparamsfiles = {
  '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/JAABA_classifier_params1.txt'
  '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/JAABA_classifier_params2.txt'
  '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/JAABA_classifier_params3.txt'
  };

CheckExperiments(expdirs,'checkprotocol',true,'checkfix',true,...
  'classifierparamsfiles',classifierparamsfiles);
