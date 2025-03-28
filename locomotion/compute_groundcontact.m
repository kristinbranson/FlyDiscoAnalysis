function [groundcontact] = compute_groundcontact(tips_velmag,varargin)

nflies = numel(tips_velmag);
ntips = size(tips_velmag{1},1);
[gc_threshold,pairs,gc_threshold_low,gc_threshold_high,minimum_bout] = myparse(varargin,'gc_threshold',20,'pairs',[],'gc_threshold_low',5,'gc_threshold_high',20,'minimum_bout',[]);

% single threshold
if numel(varargin) == 2
    for i = 1:nflies
        for j = 1:ntips
            %         currtip_velmag = tips_velmag{i}(j,:);
            currtip_velmag = tips_velmag{i}(j,:);
            idx = currtip_velmag <= gc_threshold;
            currtip_gc = zeros(1,numel(currtip_velmag));
            currtip_gc(idx) = 1;
            groundcontact{i}(j,:) = currtip_gc;
        end
    end
else
    % different thresholds for each pair of legs
    if ~isempty(pairs)
        % if using pair specific thresholding there must be a threshold for
        % each pair
        assert(size(pairs,1) == numel(gc_threshold_low),'number of low thresholds does not match number of pairs')
        assert(size(pairs,1) == numel(gc_threshold_high),'number of high thresholds does not match number of pairs')

        for i = 1:nflies
            for j = 1:ntips
                % figure out which pair
                for k = 1:size(pairs,1)
                    if intersect(j,pairs(k,:))
                        curr_gc_threshold_low = gc_threshold_low(k);
                        curr_gc_threshold_high = gc_threshold_high(k);
                        currtip_velmag = tips_velmag{i}(j,:);

                        currtip_gc = schmitt_trigger(currtip_velmag,curr_gc_threshold_low,curr_gc_threshold_high);

                        groundcontact{i}(j,:) = ~currtip_gc;
                    end
                end
            end
        end

    else
        % two thresholds
        for i = 1:nflies
            for j = 1:ntips
                %         currtip_velmag = tips_velmag{i}(j,:);
                currtip_velmag = tips_velmag{i}(j,:);

                currtip_gc = schmitt_trigger(currtip_velmag,gc_threshold_low,gc_threshold_high);

                groundcontact{i}(j,:) = ~currtip_gc;
            end
        end
    end
end

% remove bouts less than minimum bout length
if ~isempty(minimum_bout)
    for fly = 1:numel(groundcontact)
        for l = 1:6
            posts = groundcontact{fly}(l,:);
            groundcontact{fly}(l,:) = RemoveSmallBouts(posts,minimum_bout);
        end
    end
end

