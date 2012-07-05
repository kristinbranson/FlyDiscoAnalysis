%% parameters

unfixable_cats = {
  'flag_aborted_set_to_1'
  'flag_redo_set_to_1'
  'no_barcode'
  'short_video'
};

minnexpsperset_success = 2;
timestamp = '20120405';
outfilename = sprintf('RetestedLineCount_%s.tsv',timestamp);
outmatfilename = sprintf('RetestedLineCount_%s.mat',timestamp);

control_linenames = {'pBDPGAL4U','FCF_pBDPGAL4U_1500437'};

%% pull data

data = SAGEListBowlExperiments('checkflags',false,'removemissingdata',false,'screen_type','primary');

%% sort by lines, sets

[line_names,~,glineidx] = unique({data.line_name});
[set_names,~,gsetidx] = unique({data.set});
nlines = numel(line_names);
[cats,~,gcatidx] = unique({data.automated_pf_category});
ncats = numel(cats);
[~,unfixableidx] = ismember(unfixable_cats,cats);

%% 

fid = fopen(outfilename,'w');


gnsets_run = nan(1,nlines);
gnsets_success = nan(1,nlines);
gnsets_somedata = nan(1,nlines);
  
gnexps_run = nan(1,nlines);
gnexps_success = nan(1,nlines);
gnexps_somedata = nan(1,nlines);

fprintf(fid,'Line name\tN. sets run\tN. sets successful\tN. sets some data\tN. exps run\tN. exps successful\tN. exps some data\tN. exps run per set\tN. exps successful per set\tN. exps some data per set\n');

for linei = 1:nlines,
  
  if ismember(line_names{linei},control_linenames), continue; end
  
  llineidx = glineidx == linei;
  lsets = unique(gsetidx(llineidx));
  
  lnsets_run = numel(lsets);
  lnsets_success = 0;
  lnsets_somedata = 0;
  
  lnexps_run = 0;
  lnexps_success = 0;
  lnexps_somedata = 0;
  
  lnexpsperset_run = [];
  lnexpsperset_success = [];
  lnexpsperset_somedata = [];
  
  for setii = 1:lnsets_run,
    seti = lsets(setii);
    ssetidx = find(gsetidx == seti);
    snexps = numel(ssetidx);
    ismanualf = [data(ssetidx).manual_pf] == 'F';
    isautomatedf = [data(ssetidx).automated_pf] == 'F';
    scatidx = gcatidx(ssetidx);
    snfailures_total = nnz(ismanualf | isautomatedf);
    snfailures_unfixable = nnz(ismember(scatidx,unfixableidx));
    snsuccess = snexps - snfailures_total;
    snsomedata = snsuccess + snfailures_total - snfailures_unfixable;
    
    sissuccess = snsuccess >= minnexpsperset_success;
    if sissuccess,
      lnsets_success = lnsets_success + 1;
    end
    
    sissomedata = snsomedata > 0;
    if sissomedata,
      lnsets_somedata = lnsets_somedata + 1;
    end
    
    lnexps_run = lnexps_run + snexps;    
    lnexps_success = lnexps_success + snsuccess;
    lnexps_somedata = lnexps_somedata + snsomedata;

    lnexpsperset_run = [lnexpsperset_run,snexps]; %#ok<AGROW>
    lnexpsperset_success = [lnexpsperset_success,snsuccess]; %#ok<AGROW>
    lnexpsperset_somedata = [lnexpsperset_somedata,snsomedata]; %#ok<AGROW>
    
  end
  
  gnsets_run(linei) = lnsets_run;
  gnsets_success(linei) = lnsets_success;
  gnsets_somedata(linei) = lnsets_somedata;
  
  gnexps_run(linei) = lnexps_run;
  gnexps_success(linei) = lnexps_success;
  gnexps_somedata(linei) = lnexps_somedata;
  
  fprintf(fid,'%s\t%d\t%d\t%d\t%d\t%d\t%d',line_names{linei},lnsets_run,lnsets_success,lnsets_somedata,lnexps_run,lnexps_success,lnexps_somedata);
  fprintf(fid,'\t%d',lnexpsperset_run(1));
  if lnsets_run > 1,
    fprintf(fid,',%d',lnexpsperset_run(2:end));
  end
  fprintf(fid,'\t%d',lnexpsperset_success(1));
  if lnsets_run > 1,
    fprintf(fid,',%d',lnexpsperset_success(2:end));
  end
  fprintf(fid,'\t%d',lnexpsperset_somedata(1));
  if lnsets_run > 1,
    fprintf(fid,',%d',lnexpsperset_somedata(2:end));
  end
  fprintf(fid,'\n');
  
end

fclose(fid);

%%

res = struct('nsets_run',{gnsets_run},'nsets_success',{gnsets_success},'nsets_somedata',{gnsets_somedata},...
  'nexps_run',{gnexps_run},'nexps_success',{gnexps_success},'nexps_somedata',{gnexps_somedata},...
  'line_names',{line_names});

save(outmatfilename,'-struct','res');
