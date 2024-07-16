function [groundcontact] = compute_groundcontact(tips_velmag,varargin)

nflies = numel(tips_velmag);
ntips = size(tips_velmag{1},1);
[gc_threshold,gc_threshold_low,gc_threshold_high] = myparse(varargin,'gc_threshold',20,'gc_threshold_low',5,'gc_threshold_high',20);

if numel(varargin) == 2
    for i = 1:nflies
        for j = 1:ntips
            %         currtip_velmag = tips_velmag{i}(j,:);
            currtip_velmag = tips_velmag{i}(j,:);
            idx = find(currtip_velmag <= gc_threshold);
            currtip_gc = zeros(1,numel(currtip_velmag));
            currtip_gc(idx) = 1;
            groundcontact{i}(j,:) = currtip_gc;
        end
    end
else
    for i = 1:nflies
        for j = 1:ntips
            %         currtip_velmag = tips_velmag{i}(j,:);
            currtip_velmag = tips_velmag{i}(j,:);

            currtip_gc = schmitt_trigger(currtip_velmag,gc_threshold_low,gc_threshold_high);

            groundcontact{i}(j,:) = ~currtip_gc;
        end
    end
end
