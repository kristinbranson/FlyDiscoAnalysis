sizedata = SAGEGetBowlData('data_type',{'sexclassifier_diagnostics_classifier_mu_area_female','sexclassifier_diagnostics_classifier_mu_area_male'},'checkflags',false,'removemissingdata',true,'dataset','score');


%%

bad_categories = {
  'bad_number_of_flies'
  'bad_number_of_flies_per_sex'
  'flag_aborted_set_to_1'
  'flag_redo_set_to_1'
  'short_video'};
badidx = ismember({sizedata.automated_pf_category},bad_categories) | ...
  isnan([sizedata.sexclassifier_diagnostics_classifier_mu_area_female]) | ...
  isnan([sizedata.sexclassifier_diagnostics_classifier_mu_area_male]) | ...
  [sizedata.sexclassifier_diagnostics_classifier_mu_area_female] == 0 | ...
  [sizedata.sexclassifier_diagnostics_classifier_mu_area_male] == 0;

genotype = cellfun(@(x,y) [x,'__',y],{sizedata.line_name},{sizedata.effector},'UniformOutput',false);

unique_genotypes = unique(genotype(~badidx));
[~,genotypeidx] = ismember(genotype,unique_genotypes);

mean_area_female = accumarray(genotypeidx(~badidx)',...
  [sizedata(~badidx).sexclassifier_diagnostics_classifier_mu_area_female]',...
  [numel(unique_genotypes),1],...
  @mean);
std_area_female = accumarray(genotypeidx(~badidx)',...
  [sizedata(~badidx).sexclassifier_diagnostics_classifier_mu_area_female]',...
  [numel(unique_genotypes),1],...
  @(x) std(x,1));
mean_area_male = accumarray(genotypeidx(~badidx)',...
  [sizedata(~badidx).sexclassifier_diagnostics_classifier_mu_area_male]',...
  [numel(unique_genotypes),1],...
  @mean);
std_area_male = accumarray(genotypeidx(~badidx)',...
  [sizedata(~badidx).sexclassifier_diagnostics_classifier_mu_area_male]',...
  [numel(unique_genotypes),1],...
  @(x) std(x,1));
mean_area = (mean_area_female+mean_area_male)/2;
std_area = (std_area_female+std_area_male)/2;
[sortedv,order] = sort(mean_area);


clf;
plot(repmat(1:numel(unique_genotypes),[2,1]),...
  bsxfun(@plus,mean_area(order),bsxfun(@times,std_area,[-1,1]))',...
  '-','Color',[.7,.7,.7]);
hold on;
plot(1:numel(unique_genotypes),mean_area(order),'.','Color',[0,0,0]);

wtgenotypes = {'EXT_CSMH__NoEffector_0_9999'
  'EXT_CantonS_1220002__NoEffector_0_9999'
  'EXT_DL__NoEffector_0_9999'
  'FCF_attP2_1500062__NoEffector_0_9999'
  'FCF_cantons_1500002__NoEffector_0_9999'
  'pBDPGAL4U__NoEffector_0_9999'
  'pBDPGAL4U__UAS_dTrpA1_2_0002'
  'UAH_R-1220001__NoEffector_0_9999'
  'UAH_R-1220003__NoEffector_0_9999'};

idx = ismember(unique_genotypes,wtgenotypes);
[~,reorder] = sort(order);
xticks = sort(reorder(idx)');
ordered_genotypes = unique_genotypes(order);
xticklabels = ordered_genotypes(xticks);

m = regexp(xticklabels,'^GMR_(.*)_AE|D_01','tokens','once');
idx = ~cellfun(@isempty,m);
xticklabels(idx) = cellfun(@(x) x{1},m(idx),'UniformOutput',false);
m = regexp(xticklabels,'^GMR_(.*)__','tokens','once');
idx = ~cellfun(@isempty,m);
xticklabels(idx) = cellfun(@(x) x{1},m(idx),'UniformOutput',false);
m = regexp(xticklabels,'^(.*)__NoEffector','tokens','once');
idx = ~cellfun(@isempty,m);
xticklabels(idx) = cellfun(@(x) x{1},m(idx),'UniformOutput',false);
xticklabels = regexprep(xticklabels,'^EXT_','');
xticklabels = regexprep(xticklabels,'__UAS_dTrpA1_2_0002','_dTrpA1');

axisalmosttight;

set(gca,'XTick',xticks,'XTickLabel',xticklabels);
htick = rotateticklabel(gca);

textlabels = ordered_genotypes(1:5);
m = regexp(textlabels,'^GMR_(.*)_AE|D_01','tokens','once');
idx = ~cellfun(@isempty,m);
textlabels(idx) = cellfun(@(x) x{1},m(idx),'UniformOutput',false);
m = regexp(textlabels,'^GMR_(.*)__','tokens','once');
idx = ~cellfun(@isempty,m);
textlabels(idx) = cellfun(@(x) x{1},m(idx),'UniformOutput',false);
m = regexp(textlabels,'^(.*)__NoEffector','tokens','once');
idx = ~cellfun(@isempty,m);
textlabels(idx) = cellfun(@(x) x{1},m(idx),'UniformOutput',false);
textlabels = regexprep(textlabels,'^EXT_','');
textlabels = regexprep(textlabels,'__UAS_dTrpA1_2_0002','_dTrpA1');

htextl = text(1:5,sortedv(1:5),textlabels,'HorizontalAlignment','left',...
  'VerticalAlignment','middle');

textlabels = ordered_genotypes(end-4:end);
m = regexp(textlabels,'^GMR_(.*)_AE|D_01','tokens','once');
idx = ~cellfun(@isempty,m);
textlabels(idx) = cellfun(@(x) x{1},m(idx),'UniformOutput',false);
m = regexp(textlabels,'^GMR_(.*)__','tokens','once');
idx = ~cellfun(@isempty,m);
textlabels(idx) = cellfun(@(x) x{1},m(idx),'UniformOutput',false);
m = regexp(textlabels,'^(.*)__NoEffector','tokens','once');
idx = ~cellfun(@isempty,m);
textlabels(idx) = cellfun(@(x) x{1},m(idx),'UniformOutput',false);
textlabels = regexprep(textlabels,'^EXT_','');
textlabels = regexprep(textlabels,'__UAS_dTrpA1_2_0002','_dTrpA1');

htextr = text(numel(unique_genotypes)-4:numel(unique_genotypes),...
  sortedv(end-4:end),textlabels,'HorizontalAlignment','right',...
  'VerticalAlignment','middle');

%%

unique_sets = unique({sizedata(~badidx).set});
[~,setidx] = ismember({sizedata.set},unique_sets);
setmean_area_female = accumarray(setidx(~badidx)',...
  [sizedata(~badidx).sexclassifier_diagnostics_classifier_mu_area_female]',...
  [numel(unique_sets),1],...
  @mean);
setmean_area_male = accumarray(setidx(~badidx)',...
  [sizedata(~badidx).sexclassifier_diagnostics_classifier_mu_area_male]',...
  [numel(unique_sets),1],...
  @mean);
setmean_area = (setmean_area_female+setmean_area_male)/2;
set_genotype = genotype(idx1);

[unique_set_genotypes,~,setgenotypeidx] = unique(set_genotype);

figure(2);
clf;
hold on;
for i = 1:numel(unique_set_genotypes),
  
  idxcurr = find(setgenotypeidx==i);
  ncurr = numel(idxcurr);
  if ncurr <= 1,
    continue;
  end
  [idx1,idx2] = meshgrid(idxcurr,idxcurr');  
  [ism,j] = ismember(unique_genotypes{i},wtgenotypes);
  if ism,
    color = colors(j,:);
  else
    color = [0,0,0];
  end
  plot(setmean_area(idx1(:)),setmean_area(idx2(:)),'.','Color',color);
  
end

%%

fid = fopen('SortedGenotypesBySize20130717.txt','w');
for i = order',
  fprintf(fid,'%s: %f\n',unique_genotypes{i},mean_area(i));
end
fclose(fid);