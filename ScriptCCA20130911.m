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
scatter(U(:,1),U(:,2),[],zdatacluster_norm(:,stati),'.');
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


