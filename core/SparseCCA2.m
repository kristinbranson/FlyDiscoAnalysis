function [w,e,alpha,beta,mu,gamma,cor,res] = SparseCCA2(X,K,masterI,outdis,debug,sk)

%%%% CURRENT VERSION! 1.06
%
% [w,e,alpha,beta,mu,gamma,cor,res] = SCCA2(X,K,masterI,od,debug,sk)
%
% Sparse Canonical Correlation Analysis - SCCA, is a primal-dual solver for
% the CCA problem. Given primal data of a view and a dual representation of
% the second view will provide a sparse primal weight vector (for the primal
% data) and sparse feature projection (for the dual [kernel] data)
%
% Input:  X        - Primal data of view one    [m x l]
%         K        - dual data of view two      [l x l]
%         masterI  - Starting point for e       [1 x 1]
%         od       - output display true/false
%         debug    - outputs primal-dual progression true/false
%         sk       - scaling factor for mu and gamma
%
% Output: w     - sprase weight vector      [1 x m]
%         e     - sparse projct vectors     [1 x l]
%         alpha - 1'st view dual parameters [1 x m]
%         beta  - 2'nd view dual parameters [1 x l]
%         mu    - regularsation parameter   [1 x 1]
%         gamma - lagrangian (scale factor) [1 x 1]
%         cor   - correlation value         [1 x 1]
%         res   - Optimisation solution/s   [1 x 1]
%
%
% Written by David R. Hardoon 25/06/2007
% http://homepage.mac.com/davidrh/
% D.Hardoon@cs.ucl.ac.uk
%
% Further updats done - 15/04/2008
%
% Email the author any modification that are applied to the original code.
%
% No commercial usage is allowed.
%

%% checking the input
% if (nargin == 3)
%     outdis = false;
%     debug = false;
% elseif (nargin == 4)
%     debug = false;    
% elseif (nargin < 3) || (nargin > 5)
%     disp('Incorrect number of inputs');
%     help sparseSKCCA_Deflation;
%     return;
% end

%% Initialising parameters 
% Setting the infinity norm
e = zeros(size(K,2),1);
beta = e;
e(masterI) = 1;
sa_e = e;
eDifference = 0;

% So we don't need to recomput do once.
c = X*(K(:,masterI)*e(masterI));
KK = K'*K;

% More setting up initial parameters
N = size(X,1);
w = zeros(N,1);
sa_w = w;
j = ones(N,1);

% So that we do not use the e_i 
Ij = eye(size(K,2));
Ij(masterI,masterI) = 0;

% Setting initial tolerance values
etolerance = 0.01;
globtolerance = 1e-5;

% Set trade-off to half
tau = 0.5;

% Setting the mu and gamma regularsation parameters
d1 = 2*tau*(1-tau)*c; % The reminder of the equation is zero
mu = sk*mean(abs(d1));

gamma = mean(abs(2*(1-tau)^2*Ij*KK(:,masterI)*e(masterI)));

% Computing the upper bound on the w's
C = 2*mu;

% Computing inital alpha
alpha = d1 + mu*j;

% Finding alphas that break the constraints
I = find(alpha < 0 | alpha > C); 

% Selecting the violations
ta = alpha(I);
ta((ta > 0)) = ta((ta > 0)) - C;
ta = abs(ta);

[sta,stai] = sort(ta);
I = I(stai);

if length(I) > 1000
    I = I(end-999:end);
end

pI = I;

% Initial W tolerance is set
tolerance = 0.3*abs(max(alpha(I)));

if outdis
    disp(['Selected regularsation value; mu = ' num2str(mu) ', gamma = ' num2str(gamma)]);
    tic;
end

% We don't need to work on all of e
J = find(e);

% Remembering the alpha violations
preAlpha = alpha(I);

% Initially the difference will be zero
alphaDifference = abs(alpha(I) - preAlpha);

% Flag on whether to exit
exitLoop = 1;

% Loop counter
wloop = 1;

% Do we need to compute the covariance?
skipCon = 0;

% Do we need to go over all of e, to find new violations
completeE = 1;
loo = 1;

%% The loop is repeated while there are violations and we are still working
% on the alphas and that the difference between the pervious and corrent
% alphas is minimal
while((~isempty(I) && exitLoop) || (~(sum((alphaDifference > globtolerance) == 1) == 0)))
    %% set change to true so we enter the convergence on w
    change = true;   
    N = length(I);
 
    % compute the new covariance matrix if needed
    if (skipCon == 0)   
        CX = X(I,:)*X(I,:)';  
    end

    % save the previous alphas
    preAlpha = alpha(I);    

    %%% until convergence do           
    while(change)            
        %% We can exit
        change = false;

        %% Setting the update
        lefts = CX*w(I);
        
        %%% for the found alphas
        for i=1:N                       
            %% Upper and lower bounding alpha
            needtoupdate1 = false;

            if (alpha(I(i)) > C)
                alpha(I(i)) = C;                    
                needtoupdate1 = true;
            elseif (alpha(I(i)) < 0)
                alpha(I(i)) = 0;
                needtoupdate1 = true;
            else    
                %% If alpha is between the bound values
                %% shift w if needed
                if (w(I(i)) > 0)               
                    dw = (C-alpha(I(i)))/(2*tau^2*CX(i,i));
                    w(I(i)) = w(I(i)) - dw;
                elseif (w(I(i)) < 0)
                    dw = alpha(I(i))/(2*tau^2*CX(i,i));
                    w(I(i)) = w(I(i)) + dw;
                end                 
            end 

            %% Update w if need to
            if (needtoupdate1 == true)               
                %% Computing the learning rate
                learningRate = 1/(2*tau^2*CX(i,i));

               %% Updating
                firstBit = 2*tau*(1-tau)*c(I(i)) + mu - alpha(I(i));          	                
                w(I(i)) = w(I(i)) + learningRate*(firstBit - 2*tau^2*lefts(i));              
            end     

            %% Checking that w does not skip zero
            if ((sa_w(I(i)) < 0) && (w(I(i)) > 0)) || ((sa_w(I(i)) > 0) && (w(I(i)) < 0))           
                w(I(i)) = 0;
            end

            %% Computing change
            b = w(I(i))-sa_w(I(i));
            sa_w(I(i)) = w(I(i));           

            if b ~= 0
                lefts = lefts + CX(:,i)*b;          

                %% computing the new lagrangian
                alpha(I) = 2*tau*(1-tau)*c(I) + mu - 2*tau^2*lefts;
                
                %%% did we converege enough?          
                if abs(b) > tolerance
                    change = true;                     
                end
            end            
        end
      
% 		loo = loo + 1;
% 	  
%         if loo > 1000
%            % recompute alpha using the new w's and e's
%            alpha = 2*tau*(1-tau)*c + mu*j - 2*tau^2*X*(X(I,:)'*w(I)); 
%            Ii = find(alpha + globtolerance < 0 | alpha - globtolerance > C);
%            
%            if (length(Ii) ~= length(I))
%                CX = X(I,:)*X(I,:)'; 
%            end
%            
%            loo = 1;
%         end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Working now on the e's
    
	% Check whether we need to even waste time on e
	if size(K,2) > 1
        % compute all beta's (since beta are taken into account as a shadow
        % variable i.e. they are not really computed, we are able to use their
        % value as an indication of which e's are needed)
        local_beta = 2*(1-tau)^2*Ij*KK*e-2*tau*(1-tau)*Ij*K'*X(I,:)'*w(I) + gamma;

        % Find e's that need to be worked on
        J = sort([masterI; find(local_beta < 0)]);

        % Save previous e's
        preE = e(J);

        % precompute part of lagrangian update
        oneP = 2*tau*(1-tau)*Ij(J,J)'*K(:,J)'*X(I,:)'*w(I);

        % Converging over e
        change = true;
        N = size(J,1);

        % convergence over e
        while(change)          
            change = false;
            lefts = Ij(J,J)'*KK(J,J)*e(J);

            for i=1:N         
                if J(i) ~= masterI                  
                    learningRate = 1/(4*(1-tau)^2*KK(J(i),J(i)));

                    if (learningRate > 1e+3) || (learningRate < 1e-3)
                        learningRate = 1;
                    end

                    e(J(i)) = e(J(i)) + learningRate*(oneP(i) - 2*(1-tau)^2*lefts(i) +beta(J(i)) - gamma);

                    if (e(J(i)) < 0)
                        e(J(i)) = 0; 
                    elseif (e(J(i)) > 1)
                        e(J(i)) = 1;
                    else
                        beta(J(i)) = 0;
                    end

                    b = e(J(i))-sa_e(J(i));
                    sa_e(J(i)) = e(J(i));                                    

                    if b ~= 0
                        lefts = lefts + Ij(J,J)'*KK(J,J(i))*b;   

                        if abs(b) > etolerance   
                            change = true;                        
                        end
                    end
                end
            end        
        end

        % recompute c
        c = X*(K(:,J)*e(J));
        
        % check to see if there is any difference from previous (e's)
        eDifference = abs(e(J) - preE);

        % compute new tolerance values
        etolerance = 0.3*abs(max(eDifference));
        
        % bound the tolerance values
        if (etolerance == 0) || (etolerance < globtolerance)
            etolerance = globtolerance;
        end 
    end
    
    % recompute alpha using the new w's 
    alpha = 2*tau*(1-tau)*c + mu*j - 2*tau^2*X*(X(I,:)'*w(I)); 
    
    % check to see if there is any difference from previous alpha's (e's)
    alphaDifference = abs(alpha(I) - preAlpha);

    % compute new tolerance values
    tolerance = 0.3*abs(max(alphaDifference));

    if (tolerance == 0) || (tolerance < globtolerance)            
        tolerance = globtolerance;
    end

    if debug
        disp(['Loop number ' num2str(wloop)]);
        disp(['Tolerance value = ' num2str(tolerance)]);
        disp(['Error value = ' num2str(sum(alphaDifference))]);
        disp(['Etolerance value = ' num2str(etolerance)]);
        disp(['Error evalue = ' num2str(sum(eDifference))]);
    end 
    
    % Find alphas that break the constraints
    markI = I; 
    skipCon = 0;
    
    I = find(alpha + globtolerance < 0 | alpha - globtolerance > C);

    %% breakout if need to           
    if (~isempty(I)) 
        exitLoop = 1;
        %% Selecting the maximum nf violations
        ta = alpha(I);
        ta((ta > 0)) = ta((ta > 0)) - C;
        ta = abs(ta);

        %% sorting as to select the largest violations first
        [sta,stai] = sort(ta);
        I = I(stai);

        %% sanity check - do any of the violations are repeats?
        for kp=1:length(I)                     
            lc = find(I(kp) == pI);
            if lc > 0
               pI(lc) = 0;
            end
        end

        % Grab only one copy of the violations
        pI = pI((pI ~= 0));             

        %% Adding the previous I's for which w has a non zero element
        I = sort([pI((w(pI)~=0)); I]);
        if length(I) > 1000
            I = I(end-999:end);
        end

        % check to see if we need to compute the covariance matrix again
        tmp1 = sum(repmat(markI',[length(I) 1]) == repmat(I,[1 length(markI)]));
        if (sum(tmp1) == length(tmp1)) && (length(I) == length(markI))
            skipCon = 1;
        end               

        %% Saving the current index
        pI = I;        
    else
        % No violations, we can potentionally exit the algorithm
        exitLoop = 0;
        I = pI;
    end

    % update loop number
    wloop = wloop + 1;
    
%    if wloop == 5001
%        if outdis
%            disp('Need to exit, for some reason not completely converging');
%        end
%        break;
%    end
end
%%------   End of covergence algorithm
  
%%
% Compute vector length
wv = (w'*X)*(X'*w);
ev = e'*KK*e;

% normalise e
e = e/sqrt(ev);    
        
% normalise w, but check that we found something.
if sum(w ~= 0) > 0   
    w = w/sqrt(wv);     
end

% compute the optimsation error value
res = norm(tau*X'*w - (1-tau)*K*e)^2;

% compute the correlation value
cor = w'*X*K*e;     
    
if outdis         
    disp('--------------------------------------------------');

    disp(['We have ' num2str(sum(w ~= 0)) ' non zero weights']); 
    disp(['and ' num2str(sum(e ~= 0)) ' non zero dual weights'])
    disp(['Correlation = ' num2str(cor)]);
    disp(['Optimisation value = ' num2str(res)]);
    disp(['Mu = ' num2str(mu)]);
    disp(['Gamma = ' num2str(gamma)]);     
    tmp = e(masterI);
    e(masterI) = 0;
    disp(['|e|1 = ' num2str(norm(e,1))]);
    e(masterI) = tmp;
    disp(['e*KK*e = ' num2str(ev) ', w*X*X*w = ' num2str(wv)]); 
    toc; 

    disp('--------------------------------------------------');        
end
  
% For syntax clarity - end function
return
