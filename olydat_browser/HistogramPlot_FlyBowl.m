classdef HistogramPlot_FlyBowl < SingleStatTubeAveragedPlot_FlyBowl
%HistogramPlot Plot distributions of a single stat grouped by a single var
%   No interactivity. 
%

%TODO nobody likes stacked histograms. The x-ticks are always weird too.

    properties
        Name = 'Histogram';
        UsesAuxVar = false;
    end
    
    methods
        function descStr = doPlot(obj,ax,data,bstat,grp,plotCfg,~) 
            if isempty(plotCfg.grp.Name)
                % no grouping
                leg = 'All data';
            elseif iscellstr(grp)
                [grp leg] = collapsegroup(grp);
            else
                assert(isnumeric(grp))
                leg = arrayfun(@num2str,(1:max(grp))','UniformOutput',false);
            end
            nbin = max(round(numel(data)/7),4);
            tfFail = obj.detectFailure(data);
            bstat1 = bstat(~tfFail);
            grp1 = grp(~tfFail); 
            barhist_grouped(ax,bstat1,grp1,nbin,leg);
            descStr = '';
            xlabel(plotCfg.stat.PrettyName,'interpreter','none');
            ylabel('Count');
            title(sprintf('Distribution of %s',plotCfg.stat.PrettyName),'interpreter','none');
            grid(ax,'on');
        end
    end
end