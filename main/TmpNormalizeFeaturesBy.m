%% auto-select features to remove correlations for


zdatacluster_left = zdatacluster;
isselected = false(1,nstatscurr);

idxselected = nan(1,nstatscurr);
reconstructionerr = nan(1,nstatscurr);

for i = 1:nstatscurr,

  idxnotselected = find(~isselected);
  
  % zero mean
  zdatacluster_left = bsxfun(@minus,zdatacluster_left,nanmean(zdatacluster_left,1));

  % compute covariance
  tmp = zdatacluster_left;
  tmp(isnan(tmp)) = 0;
  nnotnan = double(~isnan(zdatacluster_left))'*double(~isnan(zdatacluster_left));
  cov_zdatacluster = tmp'*tmp./nnotnan;
  
  % variance
  var_zdatacluster = nanmean(zdatacluster_left.^2,1);

  err = nan(1,nstatscurr-i+1);
  for j = 1:nstatscurr-i+1,
    
    c = cov_zdatacluster(j,:) ./ var_zdatacluster(j);
    err(j) = nansum(nanmean( (bsxfun(@times,c,zdatacluster_left(:,j)) - zdatacluster_left).^2, 1),2);
    
  end
  
  
  [minerr,j] = min(err);
  
  % update zdatacluster_left: remove the reconstruction
  c = cov_zdatacluster(j,:) ./ var_zdatacluster(j);
  r = bsxfun(@times,c,zdatacluster_left(:,j));
  idxgood = ~isnan(r);
  zdatacluster_left(idxgood) = zdatacluster_left(idxgood) - r(idxgood);

  % and remove feature j
  zdatacluster_left(:,j) = [];
    
  j1 = idxnotselected(j);
  isselected(j1) = true;
  fprintf('Statistic %d: %s, err = %f\n',i,shortstatnames{j},minerr);
  
  idxselected(i) = j1;
  reconstructionerr(i) = minerr;
  
end

% uncorrelated features
zdatacluster_zeromean = bsxfun(@minus,zdatacluster,nanmean(zdatacluster,1));
zdatacluster_uncorr = nan(size(zdatacluster));
zdatacluster_uncorr(:,idxselected(1)) = zdatacluster(:,idxselected(1));
for ii = 2:nstatscurr,
  i = idxselected(ii);
  idxgood = ~isnan(zdatacluster(:,i));
  if nnz(idxgood) < ii,
    zdatacluster_uncorr(idxgood,:) = 0;
  else
    tmp = zdatacluster(idxgood,idxselected(1:ii-1));
    tmp(isnan(tmp)) = 0;
    [~,~,zdatacluster_uncorr(idxgood,i)] = regress(zdatacluster(idxgood,i),tmp);
  end
end

%%

normalize_features_by = {
  
  'fractime_flyany_frameattemptedcopulation',{'fractime_flyany_framestop'}
  'fractime_flyany_framebackup',{'fractime_flyany_framestop'}
  'fractime_flyany_framebodyturn',{'fractime_flyany_framestop'}
  'fractime_flyany_framechase',{'fractime_flyany_framestop'}
  'fractime_flyany_framecopulation',{'fractime_flyany_framestop'}
  'fractime_flyany_framecrabwalkall',{'fractime_flyany_framestop'}
  'fractime_flyany_framecrabwalkextreme',{'fractime_flyany_framestop'}
  'fractime_flyany_framejump',{'fractime_flyany_framestop'}
  'fractime_flyany_framenotanybehavior'
  'fractime_flyany_framepivotcenter',{'fractime_flyany_framestop'}
  'fractime_flyany_framepivottail',{'fractime_flyany_framestop'}
  'fractime_flyany_framerighting',{'fractime_flyany_framestop'}
  'fractime_flyany_framestop'
  'fractime_flyany_frametouch',{'fractime_flyany_framestop'}
  'fractime_flyany_framewalk',{'fractime_flyany_framestop'}
  'fractime_flyany_framewingextension',{'fractime_flyany_framestop'}
  'fractime_flyany_framewingflick',{'fractime_flyany_framestop'}
  'fractime_flyany_framewinggrooming',{'fractime_flyany_framestop'}

  'fractime_flyany_framechase_notwingextension',{'fractime_flyany_framestop','fractime_flyany_framechase','fractime_flyany_framewingextension'}
  'fractime_flyany_framestop_notwinggrooming',{'fractime_flyany_framestop','fractime_flyany_framewinggrooming'}
  'fractime_flyany_frametouch_notchase',{'fractime_flyany_frametouch','fractime_flyany_framechase'}
  'fractime_flyany_framewingextension_notchase',{'fractime_flyany_frametouch','fractime_flyany_framechase'}
  
  'fractime_flyany_framebackup_nearfly',{'fractime_flyany_framebackup','dcenter_flyany_framemove'}
  'fractime_flyany_framebackup_nearwall',{'fractime_flyany_framebackup','dist2wall_flyany_framemove'}
  'fractime_flyany_framebackup_notnearfly_notnearwall',{'fractime_flyany_framebackup','dcenter_flyany_framemove','dist2wall_flyany_framemove'}
  
  'fractime_flyany_framecrabwalkextreme_nearfly',{'fractime_flyany_framecrabwalkextreme','dcenter_flyany_framemove'}
  'fractime_flyany_framecrabwalkextreme_nearwall',{'fractime_flyany_framecrabwalkextreme','dist2wall_flyany_framemove'}
  'fractime_flyany_framecrabwalkextreme_notnearfly_notnearwall',{'fractime_flyany_framecrabwalkextreme','dcenter_flyany_framemove','dist2wall_flyany_framemove'} 
  
  'fractime_flyany_framejump_nearfly',{'fractime_flyany_framejump','dcenter_flyany_framemove'}
  'fractime_flyany_framejump_nearwall',{'fractime_flyany_framejump','dist2wall_flyany_framemove'}
  'fractime_flyany_framejump_notnearfly_notnearwall',{'fractime_flyany_framejump','dcenter_flyany_framemove','dist2wall_flyany_framemove','dist2wall_flyany_framemove'}
  
  'fractime_flyany_framemove_nearfly',{'fractime_flyany_framestop','dcenter_flyany_framemove'}
  'fractime_flyany_framemove_nearwall',{'fractime_flyany_framestop','dist2wall_flyany_framemove'}
  'fractime_flyany_framemove_notnearfly_notnearwall',{'fractime_flyany_framestop','dcenter_flyany_framemove','dist2wall_flyany_framemove'} 
    
  'fractime_flyany_framepivotcenter_nearfly',{'fractime_flyany_framepivotcenter','dcenter_flyany_framemove'}
  'fractime_flyany_framepivotcenter_nearwall',{'fractime_flyany_framepivotcenter','dist2wall_flyany_framemove'}
  'fractime_flyany_framepivotcenter_notnearfly_notnearwall',{'fractime_flyany_framepivotcenter','dcenter_flyany_framemove','dist2wall_flyany_framemove','dist2wall_flyany_framemove'}

  'fractime_flyany_framepivottail_nearfly',{'fractime_flyany_framepivottail','dcenter_flyany_framemove'}
  'fractime_flyany_framepivottail_nearwall',{'fractime_flyany_framepivottail','dist2wall_flyany_framemove'}
  'fractime_flyany_framepivottail_notnearfly_notnearwall',{'fractime_flyany_framepivottail','dcenter_flyany_framemove','dist2wall_flyany_framemove','dist2wall_flyany_framemove'}

  'fractime_flyany_framerighting_nearfly',{'fractime_flyany_framerighting','dcenter_flyany_framemove'}
  'fractime_flyany_framerighting_nearwall',{'fractime_flyany_framerighting','dist2wall_flyany_framemove'}
  'fractime_flyany_framerighting_notnearfly_notnearwall',{'fractime_flyany_framerighting','dcenter_flyany_framemove','dist2wall_flyany_framemove','dist2wall_flyany_framemove'}

  'fractime_flyany_framestop_nearfly',{'fractime_flyany_framestop','dcenter_flyany_framestop'}
  'fractime_flyany_framestop_nearwall',{'fractime_flyany_framestop','dist2wall_flyany_framestop'}
  'fractime_flyany_framestop_notnearfly_notnearwall',{'fractime_flyany_framestop','dcenter_flyany_framestop','dist2wall_flyany_framestop'}

  'fractime_flyany_framewalk_nearfly',{'fractime_flyany_framewalk','dcenter_flyany_framemove'}
  'fractime_flyany_framewalk_nearwall',{'fractime_flyany_framewalk','dist2wall_flyany_framemove'}
  'fractime_flyany_framewalk_notnearfly_notnearwall',{'fractime_flyany_framewalk','dcenter_flyany_framemove','dist2wall_flyany_framemove','dist2wall_flyany_framemove'}
  
  'fractime_flyfemale_framebackup',{'fractime_flyany_framebackup','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framebodyturn',{'fractime_flyany_framebodyturn','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framechase',{'fractime_flyany_framechase','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framecrabwalkextreme',{'fractime_flyany_framecrabwalkextreme','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framejump',{'fractime_flyany_framejump','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framenotanybehavior',{'fractime_flyany_framenotanybehavior','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framepivotcenter',{'fractime_flyany_framepivotcenter','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framepivottail',{'fractime_flyany_framepivottail','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framerighting',{'fractime_flyany_framerighting','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framestop',{'fractime_flyany_framestop'}
  'fractime_flyfemale_frametouch',{'fractime_flyany_frametouch','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framewalk',{'fractime_flyany_framewalk','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framewingflick',{'fractime_flyany_framewingflick','fractime_flyfemale_framestop'}
  'fractime_flyfemale_framewinggrooming',{'fractime_flyany_framewinggrooming','fractime_flyfemale_framestop'}
  
  'fractime_flymale_frameattemptedcopulation',{'fractime_flyany_frameattemptedcopulation','fractime_flymale_framestop'}
  'fractime_flymale_framebackup',{'fractime_flyany_framebackup','fractime_flymale_framestop'}
  'fractime_flymale_framebodyturn',{'fractime_flyany_framebodyturn','fractime_flymale_framestop'}
  'fractime_flymale_framechase',{'fractime_flyany_framechase','fractime_flymale_framestop'}
  'fractime_flymale_framecrabwalkextreme',{'fractime_flyany_framecrabwalkextreme','fractime_flymale_framestop'}
  'fractime_flymale_framejump',{'fractime_flyany_framejump','fractime_flymale_framestop'}
  'fractime_flymale_framenotanybehavior',{'fractime_flyany_framenotanybehavior','fractime_flymale_framestop'}
  'fractime_flymale_framepivotcenter',{'fractime_flyany_framepivotcenter','fractime_flymale_framestop'}
  'fractime_flymale_framepivottail',{'fractime_flyany_framepivottail','fractime_flymale_framestop'}
  'fractime_flymale_framerighting',{'fractime_flyany_framerighting','fractime_flymale_framestop'}
  'fractime_flymale_framestop',{'fractime_flyany_framestop'}
  'fractime_flymale_frametouch',{'fractime_flyany_frametouch','fractime_flymale_framestop'}
  'fractime_flymale_framewalk',{'fractime_flyany_framewalk','fractime_flymale_framestop'}
  'fractime_flymale_framewingextension',{'fractime_flyany_framewingextensions','fractime_flymale_framestop'}
  'fractime_flymale_framewingflick',{'fractime_flyany_framewingflick','fractime_flymale_framestop'}
  'fractime_flymale_framewinggrooming',{'fractime_flyany_framewinggrooming','fractime_flymale_framestop'}

  'velmag_ctr_flyany_frameany',{'fractime_flyany_framestop'}
  'velmag_ctr_flyfemale_frameany',{'velmag_ctr_flyany_frameany','fractime_flyfemale_framestop'}
  'velmag_ctr_flymale_frameany',{'velmag_ctr_flyany_frameany','fractime_flymale_framestop'}
  'velmag_ctr_flyany_framenearwall',{'velmag_ctr_flyany_frameany','fractime_flyany_framestop_nearwall'}
  'velmag_ctr_flyany_framewalk',{'velmag_ctr_flyany_frameany'}
  'velmag_ctr_flyany_framenearfly',{'velmag_ctr_flyany_frameany'}
  'velmag_ctr_flyany_framenotnearfly_notnearwall',{'velmag_ctr_flyany_frameany'}
  'velmag_ctr_flymale_framechase',{'velmag_ctr_flymale_frameany'}
  'velmag_ctr_flyany_framejump',{'velmag_ctr_flyany_frameany'}

  'du_ctr_flyany_frameany',{'velmag_ctr_flyany_frameany'}
  'du_ctr_flyany_framebackup',{'velmag_ctr_flyany_frameany'}
  
  'du_ctr_flyany_framewalk',{'velmag_ctr_flyany_framewalk'}
  'du_ctr_flyfemale_frameany',{'velmag_ctr_flyfemale_frameany'}
  'du_ctr_flymale_frameany',{'velmag_ctr_flymale_frameany'}
  'du_ctr_flymale_framechase',{'velmag_ctr_flymale_framechase'}

  'absdtheta_flyany_frameany',{'fractime_flyany_framestop'}
  'absdtheta_flyany_framepivotcenter',{'fractime_flyany_framestop','fractime_flyany_framepivotcenter'}
  'absdtheta_flyany_framepivottail',{'fractime_flyany_framestop','fractime_flyany_framepivottail'}
  'absdtheta_flyany_framebodyturn',{'fractime_flyany_framestop','fractime_flyany_framebodyturn'}
  'absdtheta_flyany_framewalk',{'fractime_flyany_framestop','fractime_flyany_framewalk'}
  
  'absdv_cor_flyany_framecrabwalkall',{'velmag_ctr_flyany_frameany','du_ctr_flyany_frameany','absdtheta_flyany_frameany'}
  'absdv_cor_flyany_framecrabwalkextreme',{'velmag_ctr_flyany_frameany','du_ctr_flyany_frameany','absdtheta_flyany_frameany'}
  'absdv_cor_flyany_framenearfly',{'velmag_ctr_flyany_framenearfly','du_ctr_flyany_frameany','absdtheta_flyany_frameany'}
  'absdv_cor_flyany_framenearwall',{'velmag_ctr_flyany_framenearwall','du_ctr_flyany_frameany','absdtheta_flyany_frameany'}
  'absdv_cor_flyany_framenotnearfly_notnearwall',{'velmag_ctr_flyany_framenotnearfly_notnearwall','du_ctr_flyany_frameany','absdtheta_flyany_frameany'}
  'absdv_cor_flyany_framewalk',{'velmag_ctr_flyany_framewalk','du_ctr_flyany_framewalk','absdtheta_flyany_frameany'}
  
  'absyaw_flyany_framemove',{'absdtheta_flyany_frameany','absdtheta_flyany_framewalk','absdv_cor_flyany_framewalk','velmag_ctr_flyany_frameany','velmag_ctr_flyany_framewalk'}
    
  'dcenter_flyany_frameany',{'fractime_flyany_framestop'}
  'dcenter_flyany_framemove',{'dcenter_flyany_frameany','fractime_flyany_framestop'}
  'dcenter_flyany_framewingflick',{'dcenter_flyany_frameany','fractime_flyany_framestop'}
  'dcenter_flyany_framestop',{'dcenter_flyany_frameany','fractime_flyany_framestop'}
  
  'dcenter_flyfemale_frameany',{'dcenter_flyany_frameany','fractime_flyfemale_framestop'}
  'dcenter_flymale_frameany',{'dcenter_flyany_frameany','fractime_flymale_framestop'}
  
  'ddcenter_flyany_framenearfly',{'dcenter_flyany_frameany','velmag_ctr_flyany_framenearfly'}
  'ddist2wall_flyany_framenearwall',{'dist2wall_flyany_frameany','velmag_ctr_flyany_framenearwall'}
  
  'ddnose2ell_flyany_framenearfly',{'dnose2ell_flyany_frameany','velmag_ctr_flyany_framenearfly'}
  'ddnose2ell_flymale_framechase',{'dnose2ell_flymale_framechase','velmag_ctr_flymale_framechase'}
    
  'dell2nose_flyany_frameany',{'dcenter_flyany_frameany','fractime_flyany_framestop'}
  'dell2nose_flymale_frameany',{'dell2nose_flyany_frameany','dcenter_flymale_frameany','fractime_flymale_framestop'}
  'dell2nose_flyfemale_frameany',{'dell2nose_flyany_frameany','dcenter_flyfemale_frameany','fractime_flyfemale_framestop'}
  'dell2nose_flyany_framewingflick',{'dell2nose_flyany_frameany','dcenter_flyany_framewingflick','fractime_flyany_framewingflick','fractime_flyany_framestop'}

  'absangle2wall_flyany_framemove',{'dist2wall_flyany_framemove','fractime_flyany_framestop','velmag_ctr_flyany_frameany'}
  'absangle2wall_flyany_framenearwall',{'dist2wall_flyany_frameany','fractime_flyany_framestop_nearwall','velmag_ctr_flyany_framenearwall'}
  'absangle2wall_flyany_framestop',{'dist2wall_flyany_framestop','fractime_flyany_framestop_nearwall'}
  
  'absanglefrom1to2_nose2ell_flyany_framenearfly',{'dnose2ell_flyany_frameany','fractime_flyany_framestop_nearfly'}
  'absanglefrom1to2_nose2ell_flyany_framestop_nearfly',{'dnose2ell_flyany_framestop','fractime_flyany_framestop_nearfly'}
  'absanglefrom1to2_nose2ell_flyany_frametouch',{'dnose2ell_flyany_frametouch','fractime_flyany_framestop_nearfly'}
  
  'absdangle2wall_flyany_framenearwall',{'dist2wall_flyany_frameany','velmag_ctr_flyany_framenearwall'}
    
  'absphidiff_nose2ell_flymale_framechase',{'dnose2ell_flymale_framechase','velmag_ctr_flymale_framechase'}
  
  'absthetadiff_nose2ell_flyany_framenearfly',{'dnose2ell_flyany_frameany','velmag_ctr_flyany_framenearfly','fractime_flyany_framechase'}
  'absthetadiff_nose2ell_flyany_frametouch',{'dnose2ell_flyany_frametouch','velmag_ctr_flyany_framenearfly','fractime_flyany_frametouch','fractime_flyany_framechase'}
  
  'absthetadiff_nose2ell_flyfemale_frametouch',{'absthetadiff_nose2ell_flyany_frametouch','dnose2ell_flyfemale_frameany','velmag_ctr_flyfemale_framenearfly','fractime_flyfemale_frametouch','fractime_flyfemale_framechase'}
  'absthetadiff_nose2ell_flymale_frametouch',{'absthetadiff_nose2ell_flyany_frametouch','dnose2ell_flymale_frameany','velmag_ctr_flymale_framenearfly','fractime_flymale_frametouch','fractime_flymale_framechase'}
  
  'absthetadiff_nose2ell_flyfemale_framenearfly',{'absthetadiff_nose2ell_flyany_framenearfly','dnose2ell_flyfemale_frameany','velmag_ctr_flyfemale_framenearfly','fractime_flyfemale_framechase'}
  'absthetadiff_nose2ell_flymale_framenearfly',{'absthetadiff_nose2ell_flyany_framenearfly','dnose2ell_flymale_frameany','velmag_ctr_flymale_framenearfly','fractime_flymale_framechase'}
  
  'absthetadiff_nose2ell_flymale_framechase',{'absthetadiff_nose2ell_flymale_framenearfly','dnose2ell_flymale_frameany','velmag_ctr_flymale_framechase','fractime_flymale_framechase'}
  'absthetadiff_nose2ell_flymale_framewingextension',{'absthetadiff_nose2ell_flymale_framenearfly','dnose2ell_flymale_frameany','velmag_ctr_flymale_framenearfly','fractime_flymale_framewingextension'}
  
  'angleonclosestfly_flyany_framenearfly',{}
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
  

  
  'veltoward_nose2ell_flyany_framenearfly'
  'veltoward_nose2ell_flymale_framechase'
  
  'wing_angle_diff_flyany_frameany'
  
  'wing_angle_diff_flyany_framewingextension'
  
  'wing_angle_imbalance_flyany_frameany'
  
  'wing_anglel_flyany_frameany'
  'wing_angler_flyany_frameany'
  'yaw_flyany_framemove'
  
  }