%% set up path


addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStats20130824.mat';
imagerydatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/ImageryData20130725.mat';
vncdatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/VNCAnnotations20130905.csv';

%% parameters

maxfraclinesmissingdata = 1;
doanatomyprocessing = false;

%% load in data

load(datafile);
if doanatomyprocessing,
  load(imagerydatafile);
end

%% compartment translations

if doanatomyprocessing,
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

%% read in the VNC annotations

fid = fopen(vncdatafile,'r');

vnc_regions = {'VNCT1','VNCT2','VNCT3','VNCAM','VNCABD','VNCSEN','VNCVEN','VNCTEC'};

vncannotations = struct('line_names',{{}});
for i = 1:numel(vnc_regions),
  vncannotations.(vnc_regions{i}) = [];
end
while true,
  
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  if isempty(s),
    continue;
  end
  ss = regexp(s,',','split');
  line = strtrim(ss{1});
  line(line=='''') = [];
  if isempty(line) || strcmp(line,'Line'),
    continue;
  end
  if numel(ss) ~= 9,
    warning('Number of elements in line >%s< ~= 9',s);
    continue;
  end
  ss(2:end) = cellfun(@(x) strtrim(x),ss(2:end),'UniformOutput',false);
  if all(cellfun(@isempty,ss(2:end))),
    continue;
  end
  isy = ismember(ss(2:end),{'Y','YY','y','yy'});
  isn = ismember(ss(2:end),{'N','NN','n','nn'});

  v = nan(1,numel(vnc_regions));
  v(isy&~isn) = 5;
  v(isn&~isy) = 0;  

  if ~strcmp(line(end-2:end),'_01'),
    line = [line,'_01'];
  end
  
  vncannotations.line_names{end+1} = line;
  for i = 1:numel(vnc_regions);
    vncannotations.(vnc_regions{i})(end+1) = v(i);
  end    
  
end

fclose(fid);

% add to linestats

[ism,idxcurr] = ismember(linestats.line_names,vncannotations.line_names);
for i = 1:numel(vnc_regions),
  fn = vnc_regions{i};
  linestats.int_manual.(fn) = nan(1,nlines);
  linestats.int_manual.(fn)(ism) = vncannotations.(fn)(idxcurr(ism));
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

%% 

statfnscurr = {
  'velmag_ctr_flyany_frameany'
  'fractime_flyany_framestop'
  'fractime_flyany_framewinggrooming'
  'fractime_flyany_framewalk'
  'fractime_flyany_framecrabwalkextreme'
  'fractime_flyany_framecrabwalkall'
  'fractime_flyany_framebackup'
  'absdtheta_flyany_frameany'
  'fractime_flyany_framepivottail'
  'fractime_flyany_framebodyturn'
  'fractime_flyany_framepivotcenter'
  'fractime_flyany_framejump'
  'fractime_flyany_framerighting'
  'dnose2ell_flyany_frameany'
  'dcenter_flyany_frameany'
  'nflies_close_flyany_frameany'
  'fractime_flyany_frametouch'
  'fractime_flyany_framechase'
  'fractime_flyany_frameattemptedcopulation'
  'wing_angle_diff_flyany_frameany'
  'fractime_flyany_framewingextension'
  'fractime_flyany_framewingflick'
  'dist2wall_flyany_frameany'
  'fractime_flyany_framenotanybehavior'
  };

%% choose some lines

line_names_curr = {...
  '.*'
  };

%% choose all lines that have some anatomy data

nmanual = zeros(1,nlines);
compartments = fieldnames(linestats.int_manual);
ncompartments = numel(compartments);
for i = 1:ncompartments,
  fn = compartments{i};
  nmanual = nmanual + ~isnan(linestats.int_manual.(fn));
end
fracmanual = nmanual / ncompartments;

lineidxcurr = fracmanual > .5;
line_names_curr = line_names(lineidxcurr);


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

anatdata = nan(nlinescurr,ncompartments);
for i = 1:ncompartments,
  anatdata(:,i) = linestats.int_manual.(compartments{i})(lineidxcurr);
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

%% hand-selected correlation removal

normalizeby = {
  'fractime_flyany_framestop'
  'velmag_ctr_flyany_frameany'
  'dcenter_flyany_frameany'
  'dist2wall_flyany_frameany'
  };

[~,idx] = ismember(normalizeby,statfnscurr);
zdatacluster_norm = nan(size(zdatacluster));
zdatacluster_norm(:,idx(1)) = zdatacluster(:,idx(1));
for ii = 2:numel(idx),

  is = idx(1:ii-1);
  j = idx(ii);
  [~,~,zdatacluster_norm(:,j)] = regress(zdatacluster(:,j),[zdatacluster(:,is),ones(nlinescurr,1)]);
  
end

X = [zdatacluster(:,idx),ones(nlinescurr,1)];
for j = setdiff(1:nstatscurr,idx),
  idxgood = ~isnan(zdatacluster(:,j));
  [~,~,zdatacluster_norm(idxgood,j)] = regress(zdatacluster(idxgood,j),X(idxgood,:));  
  fprintf('%d: %s, average abs val = %f\n',j,statfnscurr{j},mean(abs(zdatacluster_norm(idxgood,j))));
end

% fill in nans with zeros
zdatacluster_norm_nonan = zdatacluster_norm;
zdatacluster_norm_nonan(isnan(zdatacluster_norm)) = 0;

% remove some features so that the X matrix is full-rank
statidxremove_rank = ismember(statfnscurr,...
  {'max_wing_angle_flyany_framewingflick'
  'max_absdwing_angle_flyany_framewingflick'
  'duration_flyany_framewingflick'
  'dcenter_flyany_framewingflick'
  'dell2nose_flyany_framewingflick'
  'wing_anglel_flyany_frameany'
  'fractime_flyany_framechase_notwingextension'
  'wing_angle_imbalance_flyany_frameany'}');

% % here is how I selected features to make X full rank
% tmp = zdatacluster_norm_nonan;
% tmp(:,statidxremove_rank) = [];
% maxrank = rank(tmp);
% tmpnames = shortstatnames(~statidxremove_rank);
% idxcanremove = false(1,size(tmp,2));
% for i = 1:size(tmp,2),
%   [~,~,r] = regress(tmp(:,i),tmp(:,[1:i-1,i+1:size(tmp,2)]));
%   if sum(abs(r)) <= .1,
%     fprintf('%d: %s, regression residual sum = %f\n',i,tmpnames{i},sum(abs(r)));
%   end
% 
%   if rank(tmp(:,[1:i-1,i+1:size(tmp,2)])) == maxrank,
%     fprintf('%d: %s\n',i,tmpnames{i});
%     idxcanremove(i) = true;
%   end
% end


%% sparse pca

shortstatnames_pca = shortstatnames(~statidxremove_rank);
X = zdatacluster(:,~statidxremove_rank);
X(isnan(X)) = 0;
[B SD L D paths] = spca(X,[],20,0,-10,100,[],true);

for i = 1:size(B,2);
  idx = find(abs(B(:,i)) > 0);
  fprintf('\nPC %d:\n',i);
  [~,order] = sort(abs(B(idx,i)),1,'descend');
  for j = idx(order)',
    fprintf('%s: %f\n',shortstatnames_pca{j},B(j,i));
  end
end

%% for each region, find the behavior statistics that best predict it

elasticnet_alphas = [.01,.05,.25,.5,.75,.95,.99,1];
elasticnet_alphas(1) = .01;
colors = jet(nstatscurr);

B = {};
fitinfo = {};
Bbest = {};
clear fitinfobest;

for anai = 1:ncompartments,

  idxgood = ~isnan(anatdata(:,anai));
  hfig = 200+anai;
  if ishandle(hfig),
    set(0,'CurrentFigure',hfig);
  else
    figure(hfig);
  end
  clf(hfig);
  set(hfig,'Units','pixels');
  pos = get(hfig,'Position');
  pos(3:4) = [1600,500];
  set(hfig,'Position',pos);
  hax = createsubplots(2,ceil((numel(elasticnet_alphas)+1)/2),[.025,.05],hfig);

  [counts,idx] = histc(anatdata(idxgood,anai),0:1:6);
  counts(end) = [];
  weights = nan(1,nnz(idxgood));
  for j = 1:numel(counts),
    weights(idx==j) = 1/counts(j)/numel(counts);
  end
  
  
  for j = 1:numel(elasticnet_alphas),
    
    elasticnet_alpha = elasticnet_alphas(j);

%     [B{anai,j},fitinfo{anai,j}] = lasso(zdatacluster_norm_nonan(idxgood,:),anatdata(idxgood,anai),...
%       'Alpha',elasticnet_alpha,'CV',10,'PredictorNames',shortstatnames,'Standardize',false,...
%       'Weights',weights);

    i = fitinfo{anai,j}.IndexMinMSE;
    fprintf('\n%s, alpha = %f, MSE fit, error = %f:\n',compartments{anai},elasticnet_alpha,fitinfo{anai,j}.MSE(i));
    idx = find(abs(B{anai,j}(:,i)) > 0);
    [~,order] = sort(abs(B{anai,j}(idx,i)),1,'descend');
    for k = idx(order)',
      fprintf('%s: %f\n',shortstatnames{k},B{anai,j}(k,i));
    end
    i = fitinfo{anai,j}.Index1SE;
    fprintf('\n%s, alpha = %f, 1SE fit, error = %f:\n',compartments{anai},elasticnet_alpha,fitinfo{anai,j}.MSE(i));
    idx = find(abs(B{anai,j}(:,i)) > 0);
    [~,order] = sort(abs(B{anai,j}(idx,i)),1,'descend');
    for k = idx(order)',
      fprintf('%s: %f\n',shortstatnames{k},B{anai,j}(k,i));
    end
    
    %axes(hax(j));
    set(hax(j),'ColorOrder',colors);
    hold(hax(j),'on');
    firsti = fitinfo{anai,j}.IndexMinMSE;
    t = sum(abs(B{anai,j}),1);
    h = plot(hax(j),t(firsti:end),B{anai,j}(:,firsti:end),'.-');
    set(hax(j),'Color','k');
    axisalmosttight([],hax(j));
    ylim = get(hax(j),'YLim');
    xlim = get(hax(j),'XLim');
    xlim(end) = xlim(2) + .5*diff(xlim);
    set(hax(j),'XLim',xlim);
    plot(hax(j),t(fitinfo{anai,j}.IndexMinMSE)+[0,0],ylim,'w--');
    plot(hax(j),t(fitinfo{anai,j}.Index1SE)+[0,0],ylim,'w:');
    for i = 1:nstatscurr,
      text(t(firsti),B{anai,j}(i,firsti),shortstatnames{i},'Color',colors(i,:),...
        'VerticalAlignment','middle','HorizontalAlignment','left','Interpreter','none',...
        'Parent',hax(j));
    end
    title(hax(j),sprintf('%s, alpha = %f',compartments{anai},elasticnet_alpha));
   
    drawnow;
    
  end
  
  minmses = cellfun(@(x) x.MSE(x.IndexMinMSE),fitinfo(anai,:));
  onemses = cellfun(@(x) x.MSE(x.Index1SE),fitinfo(anai,:));
  bestj = argmin(minmses);
  Bbest{anai} = B{anai,bestj};
  f = fitinfo{anai,bestj};
  f.MSE_min = f.MSE(f.IndexMinMSE);
  f.MSE_1SE = f.MSE(f.Index1SE);
  [~,f.MSE_nopredictors] = weighted_mean_cov(anatdata(idxgood,anai),weights');
  fitinfobest(anai) = f;
  
  i = fitinfobest(anai).IndexMinMSE;
  prediction = zdatacluster_norm_nonan*Bbest{anai}(:,i) + fitinfobest(anai).Intercept(i);
  
  plot(hax(end),anatdata(idxgood,anai)+rand(nnz(idxgood),1)*.25,prediction(idxgood),'k.');
  axisalmosttight([],hax(end));
  axis(hax(end),'equal');
  
  drawnow;
  
end

fracgain = ([fitinfobest.MSE_nopredictors] - [fitinfobest.MSE_min])./[fitinfobest.MSE_nopredictors];

[~,order] = sort(fracgain,2,'descend');
for anai = order,
  i = fitinfobest(anai).Index1SE;
  fprintf('\n%s, alpha = %f, MSE fit, error = %f:\n',compartments{anai},fitinfobest(anai).Alpha,fitinfobest(anai).MSE(i));
  idx = find(abs(Bbest{anai}(:,i)) > 0);
  [~,order1] = sort(abs(Bbest{anai}(idx,i)),1,'descend');
  for k = idx(order1)',
    fprintf('%s: %f\n',shortstatnames{k},Bbest{anai}(k,i));
  end  
  i = fitinfobest(anai).IndexMinMSE;
  prediction = zdatacluster_norm_nonan*Bbest{anai}(:,i);
  figure(4);
  clf;
  idxgood = ~isnan(anatdata(:,anai));
  plot(anatdata(idxgood,anai)+rand(nnz(idxgood),1)*.25,prediction(idxgood),'.');
  axisalmosttight;
  axis equal;
  input(compartments{anai});
end

save('ElasticNet_PerCompartment_20130906.mat','B','fitinfo','elasticnet_alphas','compartments','statfnscurr','shortstatnames','zdatacluster_norm_nonan','anatdata');

for anai = 1:ncompartments,

  hfig = 200+anai;
  set(hfig,'Units','pixels');
  pos = get(hfig,'Position');
  pos(4) = 850;
  set(hfig,'Position',pos);
  SaveFigLotsOfWays(hfig,sprintf('ElasticNet_%s_20130906',compartments{anai}));

end


%% logistic regression


elasticnet_alphas = [1];
colors = jet(nstatscurr);

B = {};
fitinfo = {};
Bbest = {};
clear fitinfobest;

for anai = 1:ncompartments,

  idxgood = ~isnan(anatdata(:,anai)) & anatdata(:,anai) ~= 2.5;
  hfig = 200+anai;
  if ishandle(hfig),
    set(0,'CurrentFigure',hfig);
  else
    figure(hfig);
  end
  clf(hfig);
  set(hfig,'Units','pixels');
  pos = get(hfig,'Position');
  pos(3:4) = [500,850];
  set(hfig,'Position',pos);
  hax = createsubplots(2,ceil((numel(elasticnet_alphas)+1)/2),[.025,.05],hfig);

  confidence = abs(2*anatdata(idxgood,anai)/5-1);
  label = anatdata(idxgood,anai)>2.5;

  weight0 = sum(confidence(label==0));
  weight1 = sum(confidence(label==1));
  
  weights = confidence;
  weights(label==0) = weights(label==0)/weight0;
  weights(label==1) = weights(label==1)/weight1;
  weights = weights / sum(weights);  
  
  for j = 1:numel(elasticnet_alphas),
    
    elasticnet_alpha = elasticnet_alphas(j);

    [B{anai,j},fitinfo{anai,j}] = lassoglm(zdatacluster_norm_nonan(idxgood,:),label,...
      'binomial','Alpha',elasticnet_alpha,'CV',10,'PredictorNames',shortstatnames,'Standardize',false,...
      'Weights',weights,'Link','logit');

    i = fitinfo{anai,j}.IndexMinDeviance;
    fprintf('\n%s, alpha = %f, Deviance fit, error = %f:\n',compartments{anai},elasticnet_alpha,fitinfo{anai,j}.Deviance(i));
    idx = find(abs(B{anai,j}(:,i)) > 0);
    [~,order] = sort(abs(B{anai,j}(idx,i)),1,'descend');
    for k = idx(order)',
      fprintf('%s: %f\n',shortstatnames{k},B{anai,j}(k,i));
    end
    i = fitinfo{anai,j}.Index1SE;
    fprintf('\n%s, alpha = %f, 1SE fit, error = %f:\n',compartments{anai},elasticnet_alpha,fitinfo{anai,j}.Deviance(i));
    idx = find(abs(B{anai,j}(:,i)) > 0);
    [~,order] = sort(abs(B{anai,j}(idx,i)),1,'descend');
    for k = idx(order)',
      fprintf('%s: %f\n',shortstatnames{k},B{anai,j}(k,i));
    end
    
    %axes(hax(j));
    set(hax(j),'ColorOrder',colors);
    hold(hax(j),'on');
    firsti = fitinfo{anai,j}.IndexMinDeviance;
    t = sum(abs(B{anai,j}),1);
    h = plot(hax(j),t(firsti:end),B{anai,j}(:,firsti:end),'.-');
    set(hax(j),'Color','k');
    axisalmosttight([],hax(j));
    ylim = get(hax(j),'YLim');
    xlim = get(hax(j),'XLim');
    xlim(end) = xlim(2) + .5*diff(xlim);
    set(hax(j),'XLim',xlim);
    plot(hax(j),t(fitinfo{anai,j}.IndexMinDeviance)+[0,0],ylim,'w--');
    plot(hax(j),t(fitinfo{anai,j}.Index1SE)+[0,0],ylim,'w:');
    for i = 1:nstatscurr,
      text(t(firsti),B{anai,j}(i,firsti),shortstatnames{i},'Color',colors(i,:),...
        'VerticalAlignment','middle','HorizontalAlignment','left','Interpreter','none',...
        'Parent',hax(j));
    end
    title(hax(j),sprintf('%s, alpha = %f',compartments{anai},elasticnet_alpha));
   
    drawnow;
    
  end
  
  minmses = cellfun(@(x) x.Deviance(x.IndexMinDeviance),fitinfo(anai,:));
  onemses = cellfun(@(x) x.Deviance(x.Index1SE),fitinfo(anai,:));
  bestj = argmin(minmses);
  Bbest{anai} = B{anai,bestj};
  f = fitinfo{anai,bestj};
  f.Deviance_min = f.Deviance(f.IndexMinDeviance);
  f.Deviance_1SE = f.Deviance(f.Index1SE);
  fitinfobest(anai) = f;
  
  i = fitinfobest(anai).IndexMinDeviance;
  prediction = exp(zdatacluster_norm_nonan(idxgood,:)*Bbest{anai}(:,i) + fitinfobest(anai).Intercept(i));
  prediction = prediction ./ (1+prediction);
  
  centers = linspace(0,1,20);
  frac0 = hist(prediction(label==0),centers) / nnz(label==0);
  frac1 = hist(prediction(label==1),centers) / nnz(label==1);
  
  set(hax(end),'ColorOrder',[0,0,0;.7,0,0]);
  hold(hax(end),'on');
  plot(hax(end),centers,[frac0;frac1],'.-');
  axisalmosttight([],hax(end));
  
  drawnow;
  
end

for anai = 1:ncompartments,
  
  i = fitinfobest(anai).IndexMinDeviance;
  idxgood = ~isnan(anatdata(:,anai)) & anatdata(:,anai) ~= 2.5;
  confidence = abs(2*anatdata(idxgood,anai)/5-1);
  label = anatdata(idxgood,anai)>2.5;
  prediction = exp(zdatacluster_norm_nonan(idxgood,:)*Bbest{anai}(:,i) + fitinfobest(anai).Intercept(i));
  prediction = prediction ./ (1+prediction);
  err0 = sum(confidence(label==0).*prediction(label==0)) / sum(confidence(label==0));
  err1 = sum(confidence(label==1).*(1-prediction(label==1))) / sum(confidence(label==1));
  err = (err0+err1)/2;
  fitinfobest(anai).errMinDeviance = err;
  
end

[~,order] = sort([fitinfobest.errMinDeviance],2,'ascend');
for anai = order(:)',
  
  fprintf('\n%d: %s, alpha = %f, Deviance fit, error = %f:\n',anai,compartments{anai},fitinfobest(anai).Alpha,fitinfobest(anai).Deviance(i));
  idx = find(abs(Bbest{anai}(:,i)) > 0);
  [~,order1] = sort(abs(Bbest{anai}(idx,i)),1,'descend');
  for k = idx(order1)',
    fprintf('%s: %f\n',shortstatnames{k},Bbest{anai}(k,i));
  end
  
end


%% canonical correlation analysis

% remove areas without enough expression
anatidxremove_int = sum(anatdata_nonan > 2,1) <= 0;


anatdata_nonan = anatdata;
anatdata_nonan(isnan(anatdata)) = 0;
[canon_behavior,canon_anatomy,r,U,V,canon_stats] = canoncorr(zdatacluster_norm_nonan(:,~statidxremove_rank),anatdata_nonan(:,~anatidxremove_int));

tmpnames = shortstatnames(~statidxremove_rank);
fprintf('Ordered statistics in first behavior projection:\n');
[~,order] = sort(abs(canon_behavior(:,1)),1,'descend');
for i = order(1:min(numel(order),50))',
  fprintf('%s: %f\n',tmpnames{i},canon_behavior(i,1));
end

fprintf('\nOrdered regions in first anatomy projection:\n');
[~,order] = sort(abs(canon_anatomy(:,1)),1,'descend');
tmpcompartments = compartments(~anatidxremove_int);
for i = order(1:min(numel(order),50))',
  fprintf('%s: %f\n',tmpcompartments{i},canon_anatomy(i,1));
end


int_manual = linestats.int_manual;
int_manual.line_names = linestats.line_names;

hfig = 123;
figure(hfig);
clf;
stati = find(strcmp(shortstatnames,'fractime_stop'));
scatter(U(:,1),U(:,2),[],zdatacluster(:,stati),'.');
xlabel('Behavior projection 1');
ylabel('Behavior projection 2');
title('Colored by fractime stop');
hax = gca;
colorbar;
hdata_b1b2a1 = SetUpButtonDown_ReturnPointIndex(hax,U(:,1),U(:,2),{@ButtonDownFcn_ShowLineInfo,line_names_curr,'int_manual',int_manual});


hfig = 124;
figure(hfig);
clf;
h = scatter(U(:,1),U(:,2),[],V(:,1),'.');
hax = gca;
xlabel('Behavior projection 1');
ylabel('Behavior projection 2');
title('Colored by anatomy projection 1');
hdata_b1b2a1 = SetUpButtonDown_ReturnPointIndex(hax,U(:,1),U(:,2),{@ButtonDownFcn_ShowLineInfo,line_names_curr,'int_manual',int_manual});
colorbar;
%ClearButtonDown_ReturnPointIndex(hdata_b1b2a1);

hfig = 125;
figure(hfig);
clf;
plot(U(:,1),V(:,1),'k.');
xlabel('Behavior projection 1');
ylabel('Anatomy projection 1');

hfig = 126;
figure(hfig);
clf;
h = plot(V(:,1),V(:,2),'r.');
hax = gca;
xlabel('Anatomy projection 1');
ylabel('Anatomy projection 2');
axisalmosttight;
axis equal;

ax = axis;
imwidthcurr = (ax(2)-ax(1))/10;

hdata_a1bab1 = SetUpButtonDown_ReturnPointIndex(hax,V(:,1),V(:,2),...
  {@ButtonDownFcn_PlotAnatomyImageOnAxes,line_names_curr,V(:,1),V(:,2),imwidthcurr,h});

hfig = 127;
figure(hfig);
totalint = sum(anatdata_nonan(:,~anatidxremove_int),2);
clf;
h = scatter(V(:,1),V(:,2),[],totalint/ncompartments,'.');
hax = gca;
xlabel('Anatomy projection 1');
ylabel('Anatomy projection 2');
title('Colored by total intentsity');
axisalmosttight;
axis equal;

ax = axis;
imwidthcurr = 1/20;
scalecolorby = prctile(totalint,99);
set(hax,'CLim',[0,scalecolorby/ncompartments]);
colorbar;

% figure out the total weight of each compartment in each vector
anatdata_nonan_norm = bsxfun(@minus,anatdata_nonan(:,~anatidxremove_int),mean(anatdata_nonan(:,~anatidxremove_int),1));
weight1 = sum(abs(bsxfun(@times,anatdata_nonan_norm,canon_anatomy(:,1)')),1);
weight2 = sum(abs(bsxfun(@times,anatdata_nonan_norm,canon_anatomy(:,2)')),1);
[~,compartmentorder] = sort(weight1+weight2,2,'descend');
sortedcompartments = tmpcompartments(compartmentorder);

hdata_a1bab1 = SetUpButtonDown_ReturnPointIndex(hax,V(:,1),V(:,2),...
  {@ButtonDownFcn_PlotCompartmentHistogramOnAxes,line_names_curr,V(:,1),V(:,2),imwidthcurr,imwidthcurr,h,int_manual,scalecolorby,...
  'orderedcompartments',sortedcompartments});




hfig = 128;
figure(hfig);
totalint = sum(anatdata_nonan(:,~anatidxremove_int),2);
clf;
h = scatter(V(:,3),V(:,4),[],totalint/ncompartments,'.');
hax = gca;
xlabel('Anatomy projection 3');
ylabel('Anatomy projection 4');
title('Colored by total intentsity');
axisalmosttight;
axis equal;

ax = axis;
imwidthcurr = 1/20;
scalecolorby = prctile(totalint,99);
set(hax,'CLim',[0,scalecolorby/ncompartments]);
colorbar;


hdata_a1bab1 = SetUpButtonDown_ReturnPointIndex(hax,V(:,3),V(:,4),...
  {@ButtonDownFcn_PlotCompartmentHistogramOnAxes,line_names_curr,V(:,3),V(:,4),imwidthcurr,imwidthcurr,h,int_manual,scalecolorby,...
  'orderedcompartments',sortedcompartments});




%% try using cross-validation to remove features

minerr = inf;
anatidxremove_cv = false(1,ncompartments);
ncvsets = 10;
cvset = ceil(linspace(eps,ncvsets,nlinescurr));
cvset = cvset(randperm(nlinescurr));

rhocurr = nan(1,ncvsets);
for cvi = 1:ncvsets,
  [cbeh,cana] = canoncorr(zdatacluster_norm_nonan(cvset~=cvi,~statidxremove_rank),anatdata_nonan(cvset~=cvi,~anatidxremove_cv));
  pana = anatdata_nonan(cvset==cvi,~anatidxremove_cv)*cana(:,1);
  pbeh = zdatacluster_norm_nonan(cvset==cvi,~statidxremove_rank)*cbeh(:,1);
  rhocurr(cvi) = corr(pana,pbeh);
end
meanrhoprev = mean(rhocurr); 
maxmeanrhoprev = meanrhoprev;

while true,

  minmeanrho = inf;
  for compi = find(~anatidxremove_cv),

    tmpidx = anatidxremove_cv;
    tmpidx(compi) = true;
    
    rhocurr = nan(1,ncvsets);
    for cvi = 1:ncvsets,
      
      [cbeh,cana] = canoncorr(zdatacluster_norm_nonan(cvset~=cvi,~statidxremove_rank),anatdata_nonan(cvset~=cvi,~tmpidx));
      
      pana = anatdata_nonan(cvset==cvi,~tmpidx)*cana(:,1);
      pbeh = zdatacluster_norm_nonan(cvset==cvi,~statidxremove_rank)*cbeh(:,1);
      rhocurr(cvi) = corr(pana,pbeh);
      
    end
    
    meanrhocurr = mean(rhocurr);
    fprintf('%s: %f\n',compartments{compi},meanrhocurr);
    if meanrhocurr < minmeanrho,
      minmeanrho = meanrhocurr;
      worstcompi = compi;
    end
    
  end
  
  if minmeanrho < maxmeanrhoprev*.95,
    break;
  end
  
  fprintf('Removing %s, change in meanrho = %f\n',compartments{worstcompi},minmeanrho-meanrhoprev);
  meanrhoprev = minmeanrho;
  maxmeanrhoprev = max(maxmeanrhoprev,meanrhoprev);
  anatidxremove_cv(worstcompi) = true;
  
end

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