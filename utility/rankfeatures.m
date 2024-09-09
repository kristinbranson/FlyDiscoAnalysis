function [idx,z] = rankfeatures(X,group,varargin)
%RANKFEATURES ranks key features by class separability criteria
%
%   [IDX, Z] = RANKFEATURES(X,GROUP) ranks the features in X using an
%   independent evaluation criterion for binary classification, by default
%   a two-sample T-test. X is a matrix where every column is an observed
%   vector and the number of rows corresponds to the original number of
%   features. GROUP contains the class labels.  
%
%   IDX is the list of indices to the rows in X with the most significant
%   features. Z is the absolute value of the criterion used (see below). 
%
%   GROUP can be a numeric vector or a cell array of strings; numel(GROUP)
%   is the same as the number of columns in X, and numel(unique(GROUP)) is
%   equal to 2.  
%
%   RANKFEATURES(...,'CRITERION',C) sets the criterion used to assess the
%   significance of every feature for separating two labeled groups.
%   Options are
%    'ttest' (default) Absolute value two-sample T-test with pooled
%                      variance estimate
%    'entropy'         Relative entropy, also known as Kullback-Liebler
%                      distance or divergence
%    'bhattacharyya'   Minimum attainable classification error or
%                      Chernoff bound
%    'roc'             Area between the empirical receiver operating
%                      characteristic (ROC) curve and the random classifier
%                      slope
%    'wilcoxon'        Absolute value of the standarized U-statistic of a
%                      two-sample unpaired Wilcoxon test, also known as
%                      Mann-Whitney 
%    
%   Notes: 1) 'ttest', 'entropy', and 'bhattacharyya' assume normal
%   distributed classes while 'roc' and 'wilcoxon' are nonparametric tests,
%   2) all tests are feature independent.
%
%   RANKFEATURES(...,'CCWEIGHTING',ALPHA) uses correlation information to 
%   outweight the Z value of potential features using
%                         Z * (1-ALPHA*(RHO))
%   where RHO is the average of the absolute values of the cross-
%   correlation coefficient between the candidate feature and all
%   previously selected features. ALPHA sets the weighting factor. It is a
%   scalar value between 0 and 1. When ALPHA = 0 (default) potential
%   features are not weighted. A large value of RHO (close to 1) outweights
%   the significance statistic; this means that features that are highly
%   correlated with the features already picked are less likely to be
%   included in the output list.
%
%   RANKFEATURES(...,'NWEIGHTING',BETA) uses regional information to
%   outweight the Z value of potential features using
%                         Z * (1-exp(-(DIST/BETA).^2))
%   where DIST is the distance (in rows) between the candidate feature and
%   previously selected features. BETA sets the weighting factor. It is
%   greater or equal than 0. When BETA = 0 (default) potential features are
%   not weighted. A small DIST (close to 0) outweights the significance
%   statistic of only close features. This means that features that are
%   close to already picked features are less likely to be included in the
%   output list. This option is useful for extracting features from time
%   series with temporal correlation.
%
%   BETA can also be a function of the feature location, specified using @
%   or an anonymous function. In both cases RANKFEATURES passes the row
%   position of the feature to BETA() and expects back a value greater than
%   or equal to 0.  
%
%   Note: 'CCWEIGHTING' and 'NWEIGHTING' can be used together.
%
%   RANKFEATURES(...,'NUMBEROFINDICES',N) sets the number of output indices
%   in IDX. Default is the same as the number of features when ALPHA and
%   BETA are 0, or 20 otherwise.
%
%   RANKFEATURES(...,'CROSSNORM',CN) applies independent normalization
%   across the observations for every feature. Cross-normalization ensures
%   comparability among different features, although it is not always
%   necessary because some criteria already account for this.
%   Options are
%      'none' (default)  Intensities are not cross-normalized.
%      'meanvar'         x_new = (x - mean(x))/std(x)  
%      'softmax'         x_new = (1+exp((mean(x)-x)/std(x)))^-1
%      'minmax'          x_new = (x - min(x))/(max(x)-min(x))
%
%   Examples: 
%       % Find a reduced set of genes that is sufficient for
%       % differentiating breast cancer cells from all other types of
%       % cancer in the t-matrix NCI60 data set:
%       load NCI60tmatrix
%       % get a logical index vector to the breast cancer cells
%       BC = GROUP == 8;
%       % feature selection
%       I = rankfeatures(X,BC,'NumberOfIndices',12);
%       % test features with a linear discriminant classifier
%       C = classify(X(I,:)',X(I,:)',double(BC));
%       cp = classperf(BC,C);
%       cp.CorrectRate
%
%       % Use cross-correlation weighting to further reduce the required
%       % number of genes.
%       I = rankfeatures(X,BC,'CCWeighting',0.7,'NumberOfIndices',8);
%       C = classify(X(I,:)',X(I,:)',double(BC));
%       cp = classperf(BC,C);
%       cp.CorrectRate
%
%       % Find the discriminant peaks of two groups of signals with
%       % Gaussian pulses modulated by two different sources
%       load GaussianPulses
%       f = rankfeatures(y',grp,'NWeighting',@(x) x/10+5,'NumberOfIndices',5);
%       plot(t,y(grp==1,:),'b',t,y(grp==2,:),'g',t(f),1.35,'vr')
%
%   See also CANCERDETECTDEMO, CLASSIFY, CLASSPERF, CROSSVALIND,
%   RANDFEATURES, SEQUENTIALFS, SVMCLASSIFY.

%   Copyright 2003-2009 The MathWorks, Inc.


% References: 
% [1] Theodoridis, S. and Koutroumbas, K.  (1999) Pattern Recognition,
%     Academic Press, pp. 341-342.
% [2] Huan Liu, Hiroshi Motoda (1998) Feature Selection for Knowledge
%     Discovery and Data Mining, Kluwer Academic Publishers

% Example reference:
% [3] D. T. Ross, et.al. (March, 2000) Systematic Variation in Gene
%     Expression Patterns in Human Cancer Cell Lines, Published in Nature
%     Genetics, vol. 24, no. 3, pp. 227-235  


bioinfochecknargin(nargin,2,mfilename)

% set defaults
method = 'ttest';
cnorm = 'none';
alpha = 0;
beta = 0;
numIndices = NaN; % dependent, NaN will force to set it if it was not given
betaIsHandle = false;

[numPoints, numSamples] = size(X);
group = group(:);

if numel(varargin)
    if rem(numel(varargin),2)
        error(message('bioinfo:rankfeatures:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'criterion','crossnorm','ccweighting','nweighting',...
              'numberofindices'};
    for j=1:2:numel(varargin)
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:rankfeatures:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:rankfeatures:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % criterion
                    % note that we misnamed 'bhattacharyya' as
                    % 'brattacharyya' in earlier releases so keep this
                    % around for backwards compatibility.
                    methods = {'ttest','entropy','roc','wilcoxon','bhattacharyya','brattacharyya'};
                    method = strmatch(lower(pval),methods); 
                    if isempty(method) 
                      error(message('bioinfo:rankfeatures:NotValidStatistic'))
                    end
                    method = methods{method};
                case 2 % cross normalization
                    cnorms = {'none','meanvar','softmax','minmax'};
                    cnorm = strmatch(lower(pval),cnorms); 
                    if isempty(cnorm) 
                      error(message('bioinfo:rankfeatures:NotValidNormalization'))
                    end
                    cnorm = cnorms{cnorm};
                case 3 % ccweighting
                    alpha = pval(1);
                    if alpha<0 || alpha>1
                        error(message('bioinfo:rankfeatures:NotValidAlpha'))
                    end
                case 4 % nweighting
                    beta = pval;
                    betaIsHandle = isa(beta,'function_handle');
                    if ~isscalar(beta) && ~betaIsHandle
                        error('bioinfo:rankfeatures:NotValidBeta',...
                              'STEP must be a scalar or a function handler.')
                    end
                    if betaIsHandle
                        if any(beta(1:numPoints)<0)
                            error('bioinfo:rankfeatures:NotValidBeta',...
                            'Function handle BETA must return a positive number when is evaluated by any row index.')
                        end
                    else
                        if beta<0
                           error('bioinfo:rankfeatures:NotValidBeta',...
                                  'BETA must be a positive number.')
                        end
                    end
                case 5 % Number of indices
                    numIndices = round(pval(1));
                    if ~isnumeric(numIndices) || numIndices<1 || numIndices>numPoints
                        error(message('bioinfo:rankfeatures:NotValidN'))
                    end
            end
        end
    end
end
                    
% validate id and Y and some consolidation of inputs
if ~isnumeric(X) || ~isreal(X)
   error(message('bioinfo:rankfeatures:NotNumericAndReal')) 
end
if numel(group) ~= numSamples
   error(message('bioinfo:rankfeatures:NotEqualNumberOfClassLabels'))
end

group = grp2idx(group); % at this point group is numeric only, second 
                        % output of grp2idx is not needed.
todel = find(isnan(group));
if ~isempty(todel)  % remove observations not pre-classified
    X(todel,:) = [];
    numSamples = numSamples - length(todel);
end

if min(accumarray(group,1))<1 ||  max(group)~=2
    error(message('bioinfo:rankfeatures:MissingObservations'))
end

group = group==1; % change to logical
                        
% set some other dependent defaults
if isnan(numIndices)
    if alpha == 0 && ~betaIsHandle && beta==0
        numIndices = numPoints;
    else
        numIndices = 20;
    end
end

switch cnorm
    case 'meanvar'
        X=(X-repmat(nanmean(X,2),1,numSamples))./repmat(nanstd(X,[],2),1,numSamples);
    case 'softmax'
        X=1./(1+exp(-(X-repmat(nanmean(X,2),1,numSamples))./repmat(nanstd(X,[],2),1,numSamples)));
    case 'minmax'
        X=(X-repmat(min(X,[],2),1,numSamples))./repmat(max(X,[],2)-min(X,[],2),1,numSamples);
    otherwise
        % do nothing
end

n1 = sum(group);
n0 = numSamples - n1;

% do some common operations when needed
switch method
    case {'ttest','entropy','brattacharyya','bhattacharyya'}
        m1 = nanmean(X(:,group),2);
        m0 = nanmean(X(:,~group),2);
        v1 = nanvar(X(:,group),[],2);
        v0 = nanvar(X(:,~group),[],2);
end

switch method
    case 'ttest'
        % like in Zhu et.al. PNAS Dec, 2003
        z = abs((m1-m0)./ sqrt(v1/n1 + v0/n0));
    case 'entropy'
        % Theodoridis & Koutroumbas, Acad.Press, pp 152. 
        z = (v1./v0+v0./v1-2)/2+(m1-m0).^2.*(1./v1+1./v0)/2;
    case {'bhattacharyya','brattacharyya'}
        % Theodoridis & Koutroumbas, Acad.Press, pp 152. 
        s1 = sqrt(v1);s0 = sqrt(v0);
        z = ((m1-m0).^2)./(s1+s0)/4 + log((s1+s0)./sqrt(s1.*s0)/2)/2;
    case 'wilcoxon'
        % Mann-Whitney-Wilcoxon
        z = zeros(numPoints,1);
        for i = 1:numPoints
          isdata = ~isnan(X(i,:));
          if ~any(group&isdata) || ~any(~group&isdata),
            z(i) = nan;
          else
            [~,~,stats] = ranksum(X(i,group&isdata)',X(i,~group&isdata)');
            z(i) = stats.ranksum;
          end
        end
        z = abs( (z-nanmean(z)) / nanstd(z) );
        %z = abs(zscore(z));
    case 'roc'
        % Empirical Receiving Operating Characteristic Curve
        % Theodoridis & Koutroumbas, Acad.Press, pp 149.
        XallNaNs = sum(isnan(X(:,group)),2)/n1 > 0.5;
        XallNaNs = XallNaNs | (sum(isnan(X(:,~group)),2)/n0 > 0.5);
        [~,H]=sort(X,2);
        z = abs(sum(cumsum(~group(H),2).*group(H),2)/n1/n0-0.5);
        z(XallNaNs) = NaN;
    otherwise
        error(message('bioinfo:rankfeatures:NotValidStatistic'))
end

if alpha == 0 && ~betaIsHandle && beta==0
    % fast, idx to most significant markers can be found just by sorting
    [~,idx]=sort(z,1,'descend');
    % sort places NaNs above Inf, we need to move them to the back of idx
    numNaNs = sum(isnan(z));
    idx = [idx(numNaNs+1:end);idx(1:numNaNs)];
    idx = idx(1:numIndices);
else
    % this algorithm is based on the ideas presented in Theodoridis &
    % Koutroumbas, Acad.Press, pp 158-159, however the formulas are
    % different due to usability reasons, i.e. alpha and beta are easier to
    % understand and tune up.
    
    % allocate arrays and some pre-computations
    weights2 = ones(numPoints,1);
    idx = zeros(numIndices,1);
    rho = zeros(numPoints,1);
    sumXsq = sum(X.^2,2);
    [~,idx(1)]=max(z);  %most significant features
    for n = 1:numIndices-1
        % computing sums of cross-correlation coefficients
        if alpha > 0
            rho = rho + abs(sum(repmat(X(idx(n),:),numPoints,1).*X,2)./...
                            sqrt(sumXsq*sumXsq(idx(n))));
        end
        weights1 = (1-alpha*(rho/(n+1)));
        if betaIsHandle 
            % weighting out close features
            weights2 = min([weights2,...
                       1-exp(-(((1:numPoints)'-idx(n))./beta(1:numPoints)').^2)],[],2);
        elseif beta>0 % beta is just an scalar
            weights2 = min([weights2,...
                       1-exp(-(((1:numPoints)'-idx(n))./beta).^2)],[],2);
        else % only invalidate the current feature
            weights2(idx(n)) = 0;
        end
        % choose next feature
        [~,idx(n+1)] = max(z .* weights1 .* weights2); 
    end
end

