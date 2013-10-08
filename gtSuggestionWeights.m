function out = gtSuggestionWeights(jabfile,intsize)

if isstruct(jabfile),
  Q = jabfile;
else
  Q = load(jabfile,'-mat');
end
int = struct('exp',[],'flies',[],'tStart',[],'wt',[]);
out = struct('exp',[],'flies',[],'t',[],'wt',[]);

numpos = 0;
numneg = 0;
for expi = 1:numel(Q.x.gtExpDirNames),
  sfn = fullfile(Q.x.gtExpDirNames{expi},Q.x.file.scorefilename);
  load(sfn);
  for flies = 1:numel(allScores.scores)
    curt = allScores.tStart(flies):allScores.tEnd(flies);
    numpos = numpos + nnz(allScores.scores{flies}(curt)>0);
    numneg = numneg + nnz(allScores.scores{flies}(curt)<0);
  end
end
ntotal = numneg+numpos;
if numneg == 0 || numpos == 0,
  poswt = 1/ntotal;
  negwt = 1/ntotal;
else
  poswt = numneg/ntotal;
  negwt = numpos/ntotal;
end


for endx = 1:numel(Q.x.gtExpDirNames)
  sfn = fullfile(Q.x.gtExpDirNames{endx},Q.x.file.scorefilename);
  load(sfn);
  for flies = 1:numel(allScores.scores)
    curt = allScores.tStart(flies):allScores.tEnd(flies);
    if numel(curt)<intsize; continue; end
    numT = numel(curt)-intsize+1;
    int.exp(1,end+1:end+numT) = endx;
    int.flies(1,end+1:end+numT) = flies;
    int.tStart(1,end+1:end+numT) = curt(1:end-intsize+1);
    curwt = (allScores.scores{flies}(curt)<0)*negwt +(allScores.scores{flies}(curt)>0)*poswt ;
    
    cumwt = cumsum(curwt);
    sumwt = cumwt(intsize+1:end)-cumwt(1:end-intsize);
    sumwt = [cumwt(intsize) sumwt]; %#ok<AGROW>
    
    int.wt(1,end+1:end+numT) = sumwt;

    numT = numel(curt);
    out.exp(1,end+1:end+numT) = endx;
    out.flies(1,end+1:end+numT) = flies;
    out.t(1,end+1:end+numT) = curt;

    outwts = [];
    for ndx = 1:numel(curt)
      curii = (ndx-intsize+1):(ndx+intsize-1) - intsize+1;
      curii(curii<1) = [];
      curii(curii>numel(sumwt))=[];
      outwts(ndx) = sum(sumwt(curii)); %#ok<AGROW>
    end
    out.wt(1,end+1:end+numT) = outwts;
  end
end
out.wt = out.wt./sum(out.wt);

