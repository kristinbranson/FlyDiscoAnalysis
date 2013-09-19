%% set up path


addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStats20130912.mat';
imagerydatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/ImageryData20130725.mat';
vncdatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/VNCAnnotations20130905.csv';

%% load data

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

%% save results

outdatafile = sprintf('CollectedPrimaryPerFrameStatsAndAnatomy%s.mat',datestr(now,'yyyymmdd'));
save('-v7.3',outdatafile,...
  'allstats','allnframestotal','allstdstats','compartments','control_line_names','controlstd',...
  'controldateradius','exp2lineidx','expdirs','line_names','linemeans','linemetadata','linensets',...
  'linestats','main_control_line_name','metadata','ncompartments','nexps','nframestotal',...
  'nlinesmissing','nlines','nsets','nstats','set2lineidx','set_names','setidxcontrol',...
  'setiscontrol','setmeans','setmetadata','setnexps','setstats','statfns','vnc_regions','vncannotations');