classdef TornadoPlot_FlyBowl < SingleStatTubeAveragedPlot_FlyBowl
%TornadoPlot Plot distributions of a single stat grouped by a single var
%   No interactivity. The K-S statistic is used to probe for differences
%   between groups.
    
    properties
        Name = 'Tornado';
        UsesAuxVar = false;
    end
    
    methods
        function descStr = doPlot(obj,ax,data,bstat,grp,plotCfg,~) 
            if isempty(plotCfg.grp.Name)
                % no grouping
                grpname = '<All Data>';
            else
                grpname = plotCfg.grp.PrettyName;
            end
            tfFail = obj.detectFailure(data);
            bstat1 = bstat(~tfFail);
            grp1 = grp(~tfFail); %Testing out Kristin's solar powered keyboard!
            [h p grpDat] = groupdifferences(ax,bstat1,grp1,plotCfg.stat.PrettyName,grpname,'Tornado');
            % also plot the standard deviations in green
            nGrps = numel(grpDat);
            sig = nan(1,nGrps);
            mu = nan(1,nGrps);
            for i = 1:nGrps,
              mu(i) = mean(grpDat{i},1);
              sig(i) = std(grpDat{i},1,1);
            end
            holdState = get(ax,'NextPlot');
            set(ax,'NextPlot','add');
            errorbar(ax,1:nGrps,mu,sig,'r');
            set(ax,'NextPlot',holdState);
            descStr{1} = 'KS matrix:';            
            descStr{2} = 'h = ';
            descStr{3} = num2str(h);
            descStr{4} = 'p = ';
            descStr{5} = num2str(p,2);
        end
    end
end