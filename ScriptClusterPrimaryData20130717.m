%% set up path


addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStats20130614.mat';
imagerydatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/ImageryData20130725.mat';

%% parameters

maxfraclinesmissingdata = 1;

%% load in data

load(datafile);
load(imagerydatafile);

%% compartment translations

manual2autocomp = struct;
manual2autocomp.OL = {'AME','ME','LO','LOP'};
manual2autocomp.SEG = {'SOG'};
manual2autocomp.CAN = {'WED'};
manual2autocomp.BU = {'LB'};
manual2autocomp.WED = {'IVLP'};

fns = fieldnames(manual2autocomp);
fns2 = fieldnames(imdata(1).annotation.AME_L.normal);
addons = {'','_L','_R'};
for i = 1:numel(fns),
  fnm = fns{i};
  fprintf('%d: %s...\n',i,fnm);
  fna = manual2autocomp.(fnm);
  for j = 1:numel(imdata),
    imdatatmp = imdata(j);
    for k = 1:numel(addons),
      fnacurr = [fna{1},addons{k}];
      if ~isfield(imdata(j).annotation,fnacurr),
        continue;
      end
      datacurr = imdata(j).annotation.(fnacurr);
      if ~isfield(datacurr,'normal'),
        datacurr.normal = struct;
        for m = 1:numel(fns2),
          datacurr.normal.(fns2{m}) = [];
        end
      end
      if ~isfield(datacurr,'shrunk'),
        datacurr.shrunk = struct;
        for m = 1:numel(fns2),
          datacurr.shrunk.(fns2{m}) = [];
        end
      end
      imdata(j).annotation = rmfield(imdata(j).annotation,fnacurr);
      for l = 2:numel(fna),
        fnacurr = [fna{l},addons{k}];
        if ~isfield(imdata(j).annotation,fnacurr),
          continue;
        end
        for m = 1:numel(fns2),
          fn2 = fns2{m};
          if isfield(imdata(j).annotation.(fnacurr),'normal'),
            datacurr.normal.(fn2)(end+1,:) = imdata(j).annotation.(fnacurr).normal.(fn2);
          end
          if isfield(imdata(j).annotation.(fnacurr),'shrunk'),
            datacurr.shrunk.(fn2)(end+1,:) = imdata(j).annotation.(fnacurr).shrunk.(fn2);
          end
        end
        imdata(j).annotation = rmfield(imdata(j).annotation,fnacurr);
      end
      imdata(j).annotation.(fnm) = struct;
      for m = 1:numel(fns2),
        fn2 = fns2{m};
        imdata(j).annotation.(fnm).normal.(fn2) = median(datacurr.normal.(fn2),1);
        imdata(j).annotation.(fnm).shrunk.(fn2) = median(datacurr.shrunk.(fn2),1);
      end
    end
    
  end
end
  

%% combine these data

fns = fieldnames(imdata);
compartments = regexp(fns,'^(.*)_int$','tokens','once');
compartments = [compartments{:}];
ncompartments = numel(compartments);

linestats.int_manual = struct;
linestats.dist_manual = struct;
linestats.int_auto = struct;
linestats.dist_auto = struct;

for i = 1:numel(compartments),
  linestats.int_manual.(compartments{i}) = nan(1,nlines);
  linestats.dist_manual.(compartments{i}) = nan(1,nlines);
  linestats.int_auto.(compartments{i}) = nan(1,nlines);
  linestats.dist_auto.(compartments{i}) = nan(1,nlines);
end

[ism,imdata2linestats] = ismember({imdata.line},line_names);

nlinesmissing = 0;
nmissingmanual = zeros(1,nlines);
nmissingauto = zeros(1,nlines);

issymmetric = isfield(imdata(1).annotation,compartments);
islr = isfield(imdata(1).annotation,cellfun(@(x) [x,'_L'],compartments,'UniformOutput',false));

for i = 1:numel(line_names),

  idx = find(imdata2linestats == i);
  if isempty(idx),
    nmissingmanual(i) = ncompartments;
    nmissingauto(i) = ncompartments;
    continue;
  end

  
  for j = 1:numel(compartments),
    compartment = compartments{j};
    
    fnint = [compartment,'_int'];
    fndist = [compartment,'_dist'];

    intdata = [];
    distdata = [];
    for k = idx(:)',
      if isfield(imdata(k),fnint) && isfield(imdata(k),fndist),
        intdata = [intdata,imdata(k).(fnint)];
        distdata = [distdata,imdata(k).(fndist)];
      end
    end
    intdata(isnan(intdata)) = [];
    distdata(isnan(distdata)) = [];
  
    if isempty(intdata) || isempty(distdata),
      nmissingmanual(i) = nmissingmanual(i) + 1;
    else
      linestats.int_manual.(compartment)(i) = mean(intdata);
      linestats.dist_manual.(compartment)(i) = mean(distdata);
    end

    intdata = []; 
    distdata = []; 
    if issymmetric(j),

      for k = idx,
        
        if ~isfield(imdata(k).annotation,compartment),
          continue;
        end
        if isfield(imdata(k).annotation.(compartment),'shrunk'),
          if isfield(imdata(k).annotation.(compartment).shrunk,'intensity'),
            intdata = [intdata,imdata(k).annotation.(compartment).shrunk.intensity];
          end
          if isfield(imdata(k).annotation.(compartment).shrunk,'density'),          
            distdata = [distdata,imdata(k).annotation.(compartment).shrunk.density];
          end
        end
        if isfield(imdata(k).annotation.(compartment),'normal'),
          if isfield(imdata(k).annotation.(compartment).normal,'intensity'),
            intdata = [intdata,imdata(k).annotation.(compartment).normal.intensity];
          end
          if isfield(imdata(k).annotation.(compartment).normal,'density'),          
            distdata = [distdata,imdata(k).annotation.(compartment).normal.density];
          end
        end
        
      end
      
    else

      for k = idx,
        addons = {'_L','_R'};
        for l = 1:numel(addons),
          compartment_curr = [compartment,addons{l}];
          if ~isfield(imdata(k).annotation,compartment_curr),
            continue;
          end
          if isfield(imdata(k).annotation.(compartment_curr),'shrunk'),
            if isfield(imdata(k).annotation.(compartment_curr).shrunk,'intensity'),
              intdata = [intdata,imdata(k).annotation.(compartment_curr).shrunk.intensity];
            end
            if isfield(imdata(k).annotation.(compartment_curr).shrunk,'density'),
              distdata = [distdata,imdata(k).annotation.(compartment_curr).shrunk.density];
            end
          end
          if isfield(imdata(k).annotation.(compartment_curr),'normal'),
            if isfield(imdata(k).annotation.(compartment_curr).normal,'intensity'),
              intdata = [intdata,imdata(k).annotation.(compartment_curr).normal.intensity];
            end
            if isfield(imdata(k).annotation.(compartment_curr).normal,'density'),
              distdata = [distdata,imdata(k).annotation.(compartment_curr).normal.density];
            end
          end
        end
      end
  
    end
    
    intdata(isnan(intdata)) = [];
    distdata(isnan(distdata)) = [];
    
    if isempty(intdata) || isempty(distdata),
      nmissingauto(i) = nmissingauto(i) + 1;
    else
      linestats.int_auto.(compartment)(i) = mean(intdata);
      linestats.dist_auto.(compartment)(i) = mean(distdata);
    end
      
  end
end


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