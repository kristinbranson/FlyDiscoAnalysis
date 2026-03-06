function [meanim,minim,pvalue,maxpvaluesig,nlinesread] = ComputeAverageAnatomyAndPValue(line_names,varargin)

pvalue = [];
maxpvaluesig = [];
meanim = [];
minim = [];
nlinesread = 0;
defaultfilename = 'AllLinesAnatomyPValueData20130917.mat';

persistent pupperbounds pvaluetablefun highthresh lowthresh mask falsepositivedata imdata; %#ok<USENS>
[anatomydatadir,fdr,defaultfilename,imdata] = myparse(varargin,...
  'anatomydatadir','/groups/branson/bransonlab/projects/olympiad/AverageAnatomyData20140204',...
  'fdr',.1,'defaultfilename',defaultfilename,'imdata',[]);

if isempty(pvaluetablefun),
  
  if exist(defaultfilename,'file'),
    filename = defaultfilename;
  else
    filename = 'AllLinesAnatomyPValueData.mat';
  end
  [filename,pathname] = uigetfile('*.mat','Select file containing anatomy p-value data, e.g. AllLinesAnatomyPValueData20130917.mat',filename);
  if ~ischar(filename),
    return;
  end
  filename = fullfile(pathname,filename);
  load(filename,'pvaluetablefun','highthresh','lowthresh','mask','falsepositivedata','pupperbounds');
  
end

% if isempty(imdata),
%   
%   defaultfilename = 'ImageryData20130824.mat';
%   [filename,pathname] = uigetfile('*.mat','Select file containing anatomy p-value data, e.g. ImageryData20130824.mat',defaultfilename);
%   if ~ischar(filename),
%     return;
%   end
%   filename = fullfile(pathname,filename);
%   load(filename);
%   
% end

nlinesread = 0;
for i = 1:numel(line_names),
  
  fprintf('Processing data for line %d / %d...\n',i,numel(line_names));
  
  filename = fullfile(anatomydatadir,sprintf('meanim_%s.mat',line_names{i}));
    
  if ~exist(filename,'file'),
    warning('File %s does not exist.',filename);
    continue;

%     line_name = line_names{i};
%   
%     idxtmp = find(strcmp({imdata.line},line_name));
%     if isempty(idxtmp),
%       warning('File %s does not exist, and no registered images for %s exist',filename,line_name);
%       continue;
%     end
%     qi = [imdata(idxtmp).qi];
%     [minqi,j] = min(qi);
%     j = idxtmp(j);
%     if minqi > .3,
%       warning('File %s does not exist, and no registered images with qi score <= .3 for %s exist',filename,line_name);
%       continue;
%     end      
%     fprintf('%s: using stack %s, qi = %f\n',line_name,imdata(j).name,minqi);
%     datacurr = struct;
%     [datacurr.meanim] = loadRaw2Stack(imdata(j).raw_file_system_path);
%     imclass = class(datacurr.meanim);
%     if strcmpi(imclass,'uint8'),
%       scale = 2^8;
%     else
%       scale = 2^12;
%     end
%     datacurr.meanim = single(datacurr.meanim(:,:,:,2))/scale;
    
  else
    datacurr = load(filename);
  end
  if nlinesread == 0,
    [nr,nc,nz] = size(datacurr.meanim);
    countspos = zeros([nr,nc,nz]);
    countsneg = zeros([nr,nc,nz]);
    meanim = datacurr.meanim;
    minim = datacurr.meanim;
  else
    meanim = meanim + datacurr.meanim;
    minim = min(minim,datacurr.meanim);
  end
  countspos(datacurr.meanim>highthresh) = countspos(datacurr.meanim>highthresh) + 1;
  countsneg(datacurr.meanim<lowthresh) = countsneg(datacurr.meanim<lowthresh) + 1;
  nlinesread = nlinesread+1;
  
end

fprintf('Computing statistics...\n');

meanim = meanim / nlinesread;

if nargout >= 3,
  
  pvalue = ones([nr,nc,nz]);
  goodidx = (countspos+countsneg)>0;
  maskcurr = goodidx&mask;
  pvalue(maskcurr) = pvaluetablefun(countspos(maskcurr),countspos(maskcurr)+countsneg(maskcurr),pupperbounds(maskcurr));
  
end

if nargout >= 4,
  
  % nfalsepos(i) is the fraction of voxels in randomly selected lines that
  % have pvalues in [edges(i),edges(i+1))
  nsamplei = min(size(falsepositivedata.meanfrac,2),nlinesread);
  nfalsepos = falsepositivedata.meanfrac(:,nsamplei);
  counts = histc(pvalue(maskcurr),falsepositivedata.edges);
  % npos(i) is the fraction of voxels in these lines that have pvalues in
  % [edges(i),edges(i+1))
  npos = counts ./ nnz(maskcurr);
  % the expected fraction of these that are false positive is the ratio
  fdrvec = nfalsepos ./ npos;
  fdrvec(npos==0) = 0;
  fdri = find(fdrvec>=fdr,1);
  
  if isempty(fdri),
    maxpvaluesig = 1;
  else
    maxpvaluesig = falsepositivedata.edges(fdri);
  end
  
end
% 
% % Benjamini correction
% issig = falsepositivedata.edges' <= cumsum(npos)*fdr;
% if ~any(issig),
%   maxpvaluesig = 0;
% else
%   maxpvaluesig = falsepositivedata.edges(find(issig,1,'last')+1);
% end
