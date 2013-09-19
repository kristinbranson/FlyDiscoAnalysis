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