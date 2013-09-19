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

    [B{anai,j},fitinfo{anai,j}] = lasso(zdatacluster_norm_nonan(idxgood,:),anatdata(idxgood,anai),...
      'Alpha',elasticnet_alpha,'CV',10,'PredictorNames',shortstatnames,'Standardize',false,...
      'Weights',weights);

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
  %Bbest{anai} = B{anai,bestj};
  f = fitinfo{anai,bestj};
  f.Deviance_min = f.Deviance(f.IndexMinDeviance);
  f.Deviance_1SE = f.Deviance(f.Index1SE);
  %fitinfobest(anai) = f;
  
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
  
  i = fitinfobest(anai).IndexMinDeviance;
  fprintf('\n%d: %s, alpha = %f, Deviance fit, error = %f:\n',anai,compartments{anai},fitinfobest(anai).Alpha,fitinfobest(anai).errMinDeviance);
  idx = find(abs(Bbest{anai}(:,i)) > 0);
  [~,order2] = sort(abs(Bbest{anai}(idx,i)),1,'descend');
  for k = idx(order2)',
    fprintf('%s: %f\n',shortstatnames{k},Bbest{anai}(k,i));
  end
  
end

i = fitinfobest(anai).IndexMinDeviance;
idxgood = ~isnan(anatdata(:,anai)) & anatdata(:,anai) ~= 2.5;
confidence = abs(2*anatdata(idxgood,anai)/5-1);
label = anatdata(idxgood,anai)>2.5;
prediction = exp(zdatacluster_norm_nonan(idxgood,:)*Bbest{anai}(:,i) + fitinfobest(anai).Intercept(i));
prediction = prediction ./ (1+prediction);

[~,order1] = sort(prediction,2,'descend');
