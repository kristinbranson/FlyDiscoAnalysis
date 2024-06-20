function [tips_velmag] = compute_legtipvelmag(aptdata,dt,pxpermm,legtip_landmarknums)

tips_velmag = cell(1,numel(aptdata.pTrk));

nflies = numel(aptdata.pTrk);
ntips = numel(legtip_landmarknums);

for fly = 1:nflies
    for j = 1:ntips
        % landmark x,y
        x = squeeze(aptdata.pTrk{fly}(legtip_landmarknums(j),1,:))';
        y = squeeze(aptdata.pTrk{fly}(legtip_landmarknums(j),2,:))';
        % convert to mm
        x_mm = x./pxpermm;
        y_mm = y./pxpermm;
        % diff
        dx_mm = diff(x_mm);
        dy_mm = diff(y_mm);
        tips_velmag{fly}(j,:) = sqrt(dx_mm.^2+dy_mm.^2)./dt(aptdata.startframes(fly):aptdata.endframes(fly)-1)';
    end
end
