%% set up path


addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStats20130614.mat';

%% parameters

maxfraclinesmissingdata = 1;

%% load in data

load(datafile);

%% choose some statistics

statfnscurr = {
  'fractime_flyany_frameattemptedcopulation'
  'fractime_flyany_framebackup'
  'fractime_flyany_framebodyturn'
  'fractime_flyany_framechase'
  'fractime_flyany_framecopulation'
  'fractime_flyany_framecrabwalkall'
  'fractime_flyany_framecrabwalkextreme'
  'fractime_flyany_framejump'
  'fractime_flyany_framenotanybehavior'
  'fractime_flyany_framepivotcenter'
  'fractime_flyany_framepivottail'
  'fractime_flyany_framerighting'
  'fractime_flyany_framestop'
  'fractime_flyany_frametouch'
  'fractime_flyany_framewalk'
  'fractime_flyany_framewingextension'
  'fractime_flyany_framewingflick'
  'fractime_flyany_framewinggrooming'

  'fractime_flyany_framechase_notwingextension'
  'fractime_flyany_framestop_notwinggrooming'
  'fractime_flyany_frametouch_notchase'
  'fractime_flyany_framewingextension_notchase'
  
  'fractime_flyany_framebackup_nearfly'
  'fractime_flyany_framebackup_nearwall'
  'fractime_flyany_framebackup_notnearfly_notnearwall'
  
  'fractime_flyany_framecrabwalkextreme_nearfly'
  'fractime_flyany_framecrabwalkextreme_nearwall'
  'fractime_flyany_framecrabwalkextreme_notnearfly_notnearwall'
  
  
  'fractime_flyany_framejump_nearfly'
  'fractime_flyany_framejump_nearwall'
  'fractime_flyany_framejump_notnearfly_notnearwall'
  
  
  'fractime_flyany_framemove_nearfly'
  'fractime_flyany_framemove_nearwall'
  'fractime_flyany_framemove_notnearfly_notnearwall'
  
  
  'fractime_flyany_framepivotcenter_nearfly'
  'fractime_flyany_framepivotcenter_nearwall'
  'fractime_flyany_framepivotcenter_notnearfly_notnearwall'
  
  'fractime_flyany_framepivottail_nearfly'
  'fractime_flyany_framepivottail_nearwall'
  'fractime_flyany_framepivottail_notnearfly_notnearwall'
  
  
  'fractime_flyany_framerighting_nearfly'
  'fractime_flyany_framerighting_nearwall'
  'fractime_flyany_framerighting_notnearfly_notnearwall'
  
  'fractime_flyany_framestop_nearfly'
  'fractime_flyany_framestop_nearwall'
  'fractime_flyany_framestop_notnearfly_notnearwall'
  
  'fractime_flyany_framewalk_nearfly'
  'fractime_flyany_framewalk_nearwall'
  'fractime_flyany_framewalk_notnearfly_notnearwall'
  
  'fractime_flyfemale_framebackup'
  'fractime_flyfemale_framebodyturn'
  'fractime_flyfemale_framechase'
  'fractime_flyfemale_framecrabwalkextreme'
  'fractime_flyfemale_framejump'
  'fractime_flyfemale_framenotanybehavior'
  'fractime_flyfemale_framepivotcenter'
  'fractime_flyfemale_framepivottail'
  'fractime_flyfemale_framerighting'
  'fractime_flyfemale_framestop'
  'fractime_flyfemale_frametouch'
  'fractime_flyfemale_framewalk'
  'fractime_flyfemale_framewingflick'
  'fractime_flyfemale_framewinggrooming'
  
  'fractime_flymale_frameattemptedcopulation'
  'fractime_flymale_framebackup'
  'fractime_flymale_framebodyturn'
  'fractime_flymale_framechase'
  'fractime_flymale_framecrabwalkextreme'
  'fractime_flymale_framejump'
  'fractime_flymale_framenotanybehavior'
  'fractime_flymale_framepivotcenter'
  'fractime_flymale_framepivottail'
  'fractime_flymale_framerighting'
  'fractime_flymale_framestop'
  'fractime_flymale_frametouch'
  'fractime_flymale_framewalk'
  'fractime_flymale_framewingextension'
  'fractime_flymale_framewingflick'
  'fractime_flymale_framewinggrooming'
  
  'absangle2wall_flyany_framemove'
  'absangle2wall_flyany_framenearwall'
  'absangle2wall_flyany_framestop'
  
  'absanglefrom1to2_nose2ell_flyany_framenearfly'
  'absanglefrom1to2_nose2ell_flyany_framestop_nearfly'
  'absanglefrom1to2_nose2ell_flyany_frametouch'
  
  'absdangle2wall_flyany_framenearwall'
  
  'absdtheta_flyany_frameany'
  'absdtheta_flyany_framepivotcenter'
  'absdtheta_flyany_framepivottail'
  'absdtheta_flyany_framebodyturn'
  'absdtheta_flyany_framewalk'
  
  'absdv_cor_flyany_framecrabwalkall'
  'absdv_cor_flyany_framecrabwalkextreme'
  'absdv_cor_flyany_framenearfly'
  'absdv_cor_flyany_framenearwall'
  'absdv_cor_flyany_framenotnearfly_notnearwall'
  'absdv_cor_flyany_framewalk'
  
  'absphidiff_nose2ell_flymale_framechase'
  
  'absthetadiff_nose2ell_flyany_framenearfly'
  'absthetadiff_nose2ell_flyany_frametouch'
  
  'absthetadiff_nose2ell_flyfemale_frametouch'
  'absthetadiff_nose2ell_flymale_frametouch'
  
  'absthetadiff_nose2ell_flyfemale_framenearfly'
  'absthetadiff_nose2ell_flymale_framenearfly'
  
  'absthetadiff_nose2ell_flymale_framechase'
  'absthetadiff_nose2ell_flymale_framewingextension'
  
  'absyaw_flyany_framemove'
  
  'angleonclosestfly_flyany_framenearfly'
  'angleonclosestfly_flyany_framestop_nearfly'
  'angleonclosestfly_flyany_frametouch'
  
  'angleonclosestfly_flyfemale_framenearfly'
  'angleonclosestfly_flyfemale_frametouch'
  
  'angleonclosestfly_flymale_framechase'
  'angleonclosestfly_flymale_framenearfly'
  'angleonclosestfly_flymale_frametouch'
  'angleonclosestfly_flymale_framewingextension'
  
  'anglesub_flyany_frameany'
  
  'corfrac_maj_flyany_framepivotcenter'
  'corfrac_maj_flyany_framepivottail'
  
  'dangle2wall_flyany_framenearwall'
  'danglesub_flyany_framenearfly'
  
  'darea_flyany_frameany'
  
  'dcenter_flyany_frameany'
  'dcenter_flyany_framemove'
  'dcenter_flyany_framewingflick'
  'dcenter_flyany_framestop'
  
  'dcenter_flyfemale_frameany'
  'dcenter_flymale_frameany'
  
  'ddcenter_flyany_framenearfly'
  'ddist2wall_flyany_framenearwall'
  
  'ddnose2ell_flyany_framenearfly'
  'ddnose2ell_flymale_framechase'
    
  'dell2nose_flyany_frameany'
  'dell2nose_flymale_frameany'
  'dell2nose_flyfemale_frameany'
  'dell2nose_flyany_framewingflick'
  
  'dist2wall_flyany_frameany'
  'dist2wall_flyany_framemove'
  'dist2wall_flyany_framestop'
  
  'dist2wall_flyany_framewalk'
  
  'dist2wall_flyfemale_frameany'
  'dist2wall_flymale_frameany'
  
  'dmax_wing_angle_flyany_frameany'
  
  'dnose2ell_angle_30tomin30_flyany_frameany'
  'dnose2ell_angle_30tomin30_flyfemale_frameany'
  'dnose2ell_angle_30tomin30_flymale_frameany'
  
  'dnose2ell_flyany_frameany'
  
  'dnose2ell_flyany_framestop'
  'dnose2ell_flyany_frametouch'
  'dnose2ell_flyany_framewalk'
  
  'dnose2ell_flyfemale_frameany'
  'dnose2ell_flymale_frameany'
  
  'dnose2tail_flyany_frameany'
  'dnose2tail_flyany_framemove'
  'dnose2tail_flyfemale_frameany'
  'dnose2tail_flymale_frameany'
  
  'dtheta_flyany_frameany'
  
  'du_ctr_flyany_frameany'
  'du_ctr_flyany_framebackup'
  
  'du_ctr_flyany_framewalk'
  'du_ctr_flyfemale_frameany'
  'du_ctr_flymale_frameany'
  'du_ctr_flymale_framechase'
  
  'duration_flyany_framebackup'
  'duration_flyany_framebodyturn'
  'duration_flyany_framecrabwalkextreme'
  'duration_flyany_framejump'
  'duration_flyany_framemove'
  'duration_flyany_framepivotcenter'
  'duration_flyany_framepivottail'
  'duration_flyany_framerighting'
  'duration_flyany_framestop'
  'duration_flyany_frametouch'
  'duration_flyany_framewalk'
  'duration_flyany_framewingextension'
  'duration_flyany_framewingflick'
  'duration_flyany_framewinggrooming'
  
  'duration_flymale_frameattemptedcopulation'
  'duration_flymale_framechase'
  
  'dwing_angle_diff_flyany_frameany'
  
  'max_absdwing_angle_flyany_frameany'
  'max_absdwing_angle_flyany_framewingextension'
  'max_absdwing_angle_flyany_framewingflick'
  'max_absdwing_angle_flyany_framewinggrooming'
  
  'max_wing_angle_flyany_frameany'
  'max_wing_angle_flyany_framewingextension'
  'max_wing_angle_flyany_framewingflick'
  'max_wing_angle_flyany_framewinggrooming'
  
  'nflies_close_flyany_frameany'
  'nflies_close_flyany_framestop'
  'nflies_close_flyany_framewalk'
  'nflies_close_flyfemale_frameany'
  'nflies_close_flymale_frameany'
  'nflies_close_flymale_framechase'
  
  'velmag_ctr_flyany_frameany'
  'velmag_ctr_flyany_framejump'
  'velmag_ctr_flyany_framenearfly'
  'velmag_ctr_flyany_framenearwall'
  'velmag_ctr_flyany_framenotnearfly_notnearwall'
  'velmag_ctr_flyany_framewalk'
  'velmag_ctr_flyfemale_frameany'
  'velmag_ctr_flymale_frameany'
  'velmag_ctr_flymale_framechase'
  
  'veltoward_nose2ell_flyany_framenearfly'
  'veltoward_nose2ell_flymale_framechase'
  
  'wing_angle_diff_flyany_frameany'
  
  'wing_angle_diff_flyany_framewingextension'
  
  'wing_angle_imbalance_flyany_frameany'
  
  'wing_anglel_flyany_frameany'
  'wing_angler_flyany_frameany'
  'yaw_flyany_framemove'
  };

%% choose some lines

line_names_curr = {...
  '.*'
  };

%% filter these out

% stats
statidxcurr = false(numel(statfns),1);
for i = 1:numel(statfnscurr),
  statidxcurr = statidxcurr | ~cellfun(@isempty,regexp(statfns,['^',statfnscurr{i},'$'],'once'));
end
statidxcurr = find(statidxcurr);
statfnscurr = statfns(statidxcurr);
nstatscurr = numel(statidxcurr);

% lines
lineidxcurr = false(numel(line_names),1);
for i = 1:numel(line_names_curr),
  lineidxcurr = lineidxcurr | ~cellfun(@isempty,regexp(line_names,['^',line_names_curr{i},'$'],'once'))';
end
lineidxcurr = find(lineidxcurr);
line_names_curr = line_names(lineidxcurr);
nlinescurr = numel(lineidxcurr);

%% collect this data

datacluster = nan(nlinescurr,nstatscurr);
for ii = 1:nstatscurr,
  i = statidxcurr(ii);
  datacluster(:,ii) = linestats.normmeans.(statfns{i})(lineidxcurr);
end

%% z-score this data

setiscontrol = strcmp({setstats.metadata.line_name},main_control_line_name);
signorm = nan(1,nstatscurr);
for ii = 1:nstatscurr,
  i = statidxcurr(ii);
  signorm(ii) = nanstd(setstats.normmeans.(statfns{i}),1);
end

zdatacluster = bsxfun(@rdivide,datacluster,signorm);

%% remove statistics without enough data

statidxremove = find(sum(isnan(zdatacluster),1) >= nlinescurr*maxfraclinesmissingdata);
statfnscurr0 = statfnscurr;
statidxcurr0 = statidxcurr;
datacluster0 = datacluster;
zdatacluster0 = zdatacluster;

statfnscurr(statidxremove) = [];
statidxcurr(statidxremove) = [];
datacluster(:,statidxremove) = [];
zdatacluster(:,statidxremove) = [];
signorm(statidxremove) = [];
nstatscurr = numel(statidxcurr);

%% create short names for stats and lines for plotting

shortstatnames = statfnscurr;
shortstatnames = regexprep(shortstatnames,'_flyany','');
shortstatnames = regexprep(shortstatnames,'^(.*)_fly(.*)_(.*)','$1_$3_$2');
shortstatnames = regexprep(shortstatnames,'^fractime_frame','fractime_');
shortstatnames = regexprep(shortstatnames,'^duration_frame','duration_');
shortstatnames = regexprep(shortstatnames,'_frameany','');
shortstatnames = regexprep(shortstatnames,'frame','');

shortlinenames = line_names_curr;
shortlinenames = regexprep(shortlinenames,'GMR_','R');
shortlinenames = regexprep(shortlinenames,'_AE_01','');
shortlinenames = regexprep(shortlinenames,'_AD_01','D');

%% compute pairwise distance between lines, ignoring entries for which either has nan

% L1 distance
lined = zeros(nlinescurr,nlinescurr);
for linei = 1:nlinescurr,
  for linej = linei+1:nlinescurr,
    dcurr = abs(zdatacluster(linei,:)-zdatacluster(linej,:));
    lined(linei,linej) = nanmean(dcurr);
    lined(linej,linei) = lined(linei,linej);
  end
end
linedvec = squareform(lined,'tovector');

%% compute pairwise distance between stats, ignoring entries for which either has nan

% 1- abs(correlation coeff)
statd = zeros(nstatscurr,nstatscurr);
for stati = 1:nstatscurr,
  ignorei = isnan(zdatacluster(:,stati));
  for statj = stati+1:nstatscurr,
    ignorecurr = ignorei | isnan(zdatacluster(:,statj));
    if all(ignorecurr),
      statd(stati,statj) = nan;
    else
      r = corrcoef(zdatacluster(~ignorecurr,stati),zdatacluster(~ignorecurr,statj));
      statd(stati,statj) = 1 - abs(r(1,2));
    end
    statd(statj,stati) = statd(stati,statj);
  end
end
statdvec = squareform(statd,'tovector');

%% my version of a clustergram

cgobj = clustergram(zdatacluster',...
  'RowLabels',shortstatnames,...
  'ColumnLabels',shortlinenames,...
  'Standardize','none',...
  'Cluster','all',...
  'RowPDist',statdvec,...
  'ColumnPDist',linedvec,...
  'Linkage','average',...
  'OptimalLeafOrder',true,...
  'ImputeFun',@ClustergramImputeFun);