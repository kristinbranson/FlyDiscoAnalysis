classdef CorrelationPlot_FlyBowl < OlyDat.Plot
%CorrelationPlot Plot correlation between a single stat and an aux stat
%   Color using a single grouping var. Experiment selection using
%   XYPlotClickHandler (single-click, or bounding-box selection, with or
%   without shift-click). A correlation coefficient/p-value is computed.

    properties
        Name = 'Correlation';
        UsesAuxVar = true;
        AllowsMultiStats = false;
    end
    
    methods
        function str = doPlot(obj,ax,data,bstat,grp,pltCfg,expCoordinator) %#ok<MANU>
            assert(isscalar(pltCfg.aux));
            aux = [data.(pltCfg.aux.DataFieldName)]';
            eid = [data.experiment_id]';

            % convert string grp to numeric
            if isempty(pltCfg.grp.Name)
                % no grouping
                grpleg = {'All data'};
            else
                [grp grpleg] = collapsegroup(grp);
                if isnumeric(grpleg)
                    grpleg = arrayfun(@num2str,grpleg,'UniformOutput',false);
                end
            end
            
            reset(ax);
            zlclCorrPlot(ax,aux,bstat,grp,pltCfg.aux.PrettyName,pltCfg.stat.PrettyName,grpleg);

            assert(isequal(size(aux),size(bstat),size(eid)));
            OlyDat.XYPlotClickHandler(ax,aux,bstat,eid,expCoordinator,true);
            str = '';
        end
        
        function data = preprocessData(obj,data,pltCfg)
            ystat = pltCfg.stat;
            assert(isscalar(ystat));
            xstat = pltCfg.aux;
            assert(isscalar(xstat));
            
            stats = [xstat;ystat];
            tfIsScore = [stats.IsScore]';
            if any(tfIsScore)
                % at least one bstat is a score; need to flatten the scores
                % to the top-level
                scoreStats = stats(tfIsScore);
                assert(~any(isfield(data,{scoreStats.DataFieldName}')));
                flatdata = flattendata(data,[],...
                    @(scoreArr)obj.extractScores(scoreArr,scoreStats));
            else
                % both x, y stats are top-level
                flatdata = data;
            end
            
            assert(isequal(size(flatdata),size(data)));
            data = flatdata;            
        end
        
        % all stats should be scores. 
        function s = extractScores(obj,scoreArr,stats) %#ok<MANU>
            s = struct();
            for c = 1:numel(stats)
                st = stats(c);
                assert(st.IsScore);
                v = nan(6,1);
                for d = 1:6
                    ss = scoreArr{st.TempIdx,st.Sequence,d};
                    if ~isempty(ss)
                        v(d) = ss.(st.Name);
                    end
                end
                s.(st.DataFieldName) = nanmean(v);
            end
        end
        
    end
end

function [r p] = zlclCorrPlot(ax,x,y,g,xVar,yVar,grpleg)
switch xVar
  case {'experiment date time','exp_datenum','edt','edt_day','edt_week'}
    xtickLblFcn = 'yymmmdd';
  case {'edt_month'},
    xtickLblFcn = 'yymmm';
  otherwise
    xtickLblFcn = [];
end
[r p] = scattercorrelate(ax,x,y,g,grpleg,'',xVar,yVar,xtickLblFcn);
grid(ax,'on');
title(ax,sprintf('%s vs %s: r=%.3f, p=%.3f',yVar,xVar,r(1,2),p(1,2)),'interpreter','none');
end
