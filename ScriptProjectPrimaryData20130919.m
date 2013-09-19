%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStatsAndAnatomy20130912.mat';

%% load in data

load(datafile);

%% parameters

maxfraclinesmissingdata = 1;
doanatomyprocessing = false;
statfnset = 'many';
lineset = 'all';
disttransform = 'linearthenlog';
linear2loginflectionpt = 3;

%% choose some statistics

switch statfnset,
  
  case 'many',

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
    
  case 'few'

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
    
end

%% choose some lines

switch lineset,
  
  case 'all',

    line_names_curr = {...
      '.*'
      };
    
  case 'hasanatomyannotation',

    % choose all lines that have some anatomy data
    
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

end
    

%% filter these out

% stats
statidxcurr = false(1,numel(statfns));
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

ncompartments = numel(compartments);
anatdata = nan(nlinescurr,ncompartments);
for i = 1:ncompartments,
  anatdata(:,i) = linestats.int_manual.(compartments{i})(lineidxcurr);
end

%% hand-selected correlation removal

normalizeby = {
  'fractime_flyany_framestop'
  'velmag_ctr_flyany_frameany'
  'dcenter_flyany_frameany'
  'dist2wall_flyany_frameany'
  };

setiscontrol = strcmp({setstats.metadata.line_name},main_control_line_name);
setdata = nan(nnz(setiscontrol),nstatscurr);
for ii = 1:nstatscurr,
  i = statidxcurr(ii);
  setdata(:,ii) = setstats.normmeans.(statfns{i})(setiscontrol);
end

if isempty(normalizeby),
  datacluster_norm = datacluster;
  setdata_norm = setdata;
else
  [~,idx] = ismember(normalizeby,statfnscurr);
  datacluster_norm = nan(size(datacluster));
  datacluster_norm(:,idx(1)) = datacluster(:,idx(1));
  j = idx(1);
  fprintf('%d: %s, average abs val = %f\n',j,statfnscurr{j},mean(abs(datacluster_norm(:,j))));

  normcoeffs = cell(1,nstatscurr);
  for ii = 2:numel(idx),
    
    is = idx(1:ii-1);
    j = idx(ii);
    [normcoeffs{j},~,datacluster_norm(:,j)] = regress(datacluster(:,j),[datacluster(:,is),ones(nlinescurr,1)]);
    fprintf('%d: %s, average abs val = %f\n',j,statfnscurr{j},mean(abs(datacluster_norm(:,j))));
    
  end

  X = [datacluster(:,idx),ones(nlinescurr,1)];
  for j = setdiff(1:nstatscurr,idx),
    idxgood = ~isnan(datacluster(:,j));
    [normcoeffs{j},~,datacluster_norm(idxgood,j)] = regress(datacluster(idxgood,j),X(idxgood,:));
    fprintf('%d: %s, average abs val = %f\n',j,statfnscurr{j},mean(abs(datacluster_norm(idxgood,j))));
  end
  
  % also normalize the control set stats so that we can z-score

  setdata_norm = nan(size(setdata));

  % first feature
  j = idx(1);
  setdata_norm(:,j) = setdata(:,j);

  % next features
  for ii = 2:numel(idx),
    is = idx(1:ii-1);
    j = idx(ii);
    pred = [setdata(:,is),ones(nnz(setiscontrol),1)]*normcoeffs{j};
    setdata_norm(:,j) = setdata(:,j) - pred;
  end

  % rest of features
  X = [setdata(:,idx),ones(nnz(setiscontrol),1)];
  for j = setdiff(1:nstatscurr,idx),
    pred = X*normcoeffs{j};
    setdata_norm(:,j) = setdata(:,j) - pred;
  end

end
  
% compute mean and standard deviations
munorm = nanmean(setdata_norm,1);
signorm = nanstd(setdata_norm,1,1);

% z-score the line data
zdatacluster_norm = bsxfun(@rdivide,bsxfun(@minus,datacluster_norm,munorm),signorm);

% and the control set data
zcontrol_setdata_norm = bsxfun(@rdivide,bsxfun(@minus,setdata_norm,munorm),signorm);

% % remove nans
% zdatacluster_norm_nonan = zdatacluster_norm;
% zdatacluster_norm_nonan(isnan(zdatacluster_norm)) = 0;
% 
% % same for the control set data
% zcontrol_setdata_norm_nonan = zcontrol_setdata_norm;
% zcontrol_setdata_norm_nonan(isnan(zcontrol_setdata_norm)) = 0;

% 
% % fill in nans with zeros
% zdatacluster_norm_nonan = zdatacluster_norm;
% zdatacluster_norm_nonan(isnan(zdatacluster_norm)) = 0;
% 
% % remove some features so that the X matrix is full-rank
% statidxremove_rank = ismember(statfnscurr,...
%   {'max_wing_angle_flyany_framewingflick'
%   'max_absdwing_angle_flyany_framewingflick'
%   'duration_flyany_framewingflick'
%   'dcenter_flyany_framewingflick'
%   'dell2nose_flyany_framewingflick'
%   'wing_anglel_flyany_frameany'
%   'fractime_flyany_framechase_notwingextension'
%   'wing_angle_imbalance_flyany_frameany'}');

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

%% remove statistics without enough data

statidxremove = find(sum(isnan(zdatacluster_norm),1) >= nlinescurr*maxfraclinesmissingdata);
statfnscurr0 = statfnscurr;
statidxcurr0 = statidxcurr;
datacluster0 = datacluster;
zdatacluster_norm0 = zdatacluster_norm;
zcontrol_setdata_norm0 = zcontrol_setdata_norm;

statfnscurr(statidxremove) = [];
statidxcurr(statidxremove) = [];
datacluster(:,statidxremove) = [];
zdatacluster_norm(:,statidxremove) = [];
zcontrol_setdata_norm(:,statidxremove) = [];
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

%% non-linear transformation

if strcmpi(disttransform,'linearthenlog'),
  
  % line data
  zdatacluster_transform = nan(size(zdatacluster_norm));
  idx = abs(zdatacluster_norm)<=linear2loginflectionpt;
  zdatacluster_transform(idx) = zdatacluster_norm(idx);
  zdatacluster_transform(~idx) = (linear2loginflectionpt+...
    log(abs(zdatacluster_norm(~idx))-linear2loginflectionpt+1)).*...
    sign(zdatacluster_norm(~idx));

  % control set data
  zcontrol_setdata_transform = nan(size(zcontrol_setdata_norm));
  idx = abs(zcontrol_setdata_norm)<=linear2loginflectionpt;
  zcontrol_setdata_transform(idx) = zcontrol_setdata_norm(idx);
  zcontrol_setdata_transform(~idx) = (linear2loginflectionpt+...
    log(abs(zcontrol_setdata_norm(~idx))-linear2loginflectionpt+1)).*...
    sign(zcontrol_setdata_norm(~idx));
  
elseif strcmpi(disttransform,'log'),

  zdatacluster_transform = log(abs(zdatacluster_norm)).*sign(zdatacluster_norm);

  zcontrol_setdata_transform = log(abs(zcontrol_setdata_norm)).*sign(zcontrol_setdata_norm);

else
  
  zdatacluster_transform = zdatacluster_norm;

  zcontrol_setdata_transform = zcontrol_setdata_norm;

end

%% PCA

X = zdatacluster_transform;
X(isnan(X)) = 0;
mu_pca = mean(X,1);
X = bsxfun(@minus,X,mu_pca);
[coeff,proj_pca,latent,tsquared] = princomp(X);

% find the projection of all the control sets
Xcontrol = zcontrol_setdata_transform;
Xcontrol(isnan(Xcontrol)) = 0;
Xcontrol = bsxfun(@minus,Xcontrol,mu_pca);
projcontrol_pca = Xcontrol*coeff;

% control variance
Scontrol_pca = cov(projcontrol_pca,1);
mucontrol_pca = mean(projcontrol_pca,1);
[a,b,theta] = cov2ell(Scontrol_pca(1:2,1:2));

hfig = 1;
figure(hfig);
clf;

maxnstdsplot = 5;
colors = gray(maxnstdsplot);
for i = maxnstdsplot:-1:1,
  hcov(i) = ellipsedrawpatch(a*i/2,b*i/2,mucontrol_pca(1),mucontrol_pca(2),theta,colors(i,:));
  hold on;
end

plot(projcontrol_pca(:,1),projcontrol_pca(:,2),'kx');
hold on;
plot(proj_pca(:,1),proj_pca(:,2),'.','Color',[.8,0,0]);

%% sparse pca

shortstatnames_pca = shortstatnames(~statidxremove_rank);
X = zdatacluster_norm(:,~statidxremove_rank);
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