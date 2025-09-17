function [CoM_stabilty] = compute_CoMstability(aptdata,legtip_landmarknums,groundcontact)

debug = false;
pTrk = aptdata.pTrk;



parfor fly = 1:numel(pTrk)
    % compute CoM
    leftshoulder = squeeze(pTrk{fly}(4,:,:));
    rightshoulder  = squeeze(pTrk{fly}(5,:,:));
    ctr = (leftshoulder + rightshoulder) /2;
    % head = squeeze(pTrk{fly}(1,:,:));
    % tail = squeeze(pTrk{fly}(7,:,:));
    notum = squeeze(pTrk{fly}(6,:,:));

    % CoM_body{fly} = head + (tail-head).*0.4;
    CoM_thorax{fly} = ctr + (notum - ctr).*0.6;

    % compute CoM Stability based on doi.org/10.7554/eLife.65878
    nframes = size(CoM_thorax{fly}, 2);
    CoM_stabilty{fly} = nan(1, nframes);

    % find frames with 3 or more feet on the ground
    nfeet = sum(groundcontact{fly},1);
    valid_frames = find(nfeet >= 3);

    if isempty(valid_frames)
        continue;
    end

    for i = 1:numel(valid_frames)
        frame_idx = valid_frames(i);

        % find legs with groun contact
        idx_gc = groundcontact{fly}(:,frame_idx);
        legs_idx = legtip_landmarknums(idx_gc);
        % stance legs polygon vertices
        XV = pTrk{fly}(legs_idx,1,frame_idx);
        YV = pTrk{fly}(legs_idx,2,frame_idx);
        % CoM
        X_CoM = CoM_thorax{fly}(1,frame_idx);
        Y_CoM = CoM_thorax{fly}(2,frame_idx);
        % outside polygon
        if ~inpolygon(X_CoM,Y_CoM,XV,YV)
            CoM_stabilty{fly}(frame_idx) = 0;
        else
            % find the shortest distance to the polygon when inside
            % claude.ai code
            CoM_stabilty{fly}(frame_idx) = fastPointToPolygonDistance([X_CoM, Y_CoM], [XV, YV]);
        end
    end


end

if debug
    h = figure
    hold on
    fly = 1;
    frm = 1000;
    plot(pTrk{fly}(:,1,frm),pTrk{fly}(:,2,frm),'.')
    plot(ctr(1,frm),ctr(2,frm),'.r')
    plot(CoM_body{fly}(1,frm),CoM_body{fly}(2,frm),'ok')
    plot(CoM_thorax{fly}(1,frm),CoM_thorax{fly}(2,frm),'*r')
end


end

function minDist = fastPointToPolygonDistance(point, polygon)
% Optimized version using vectorized operations

% Close polygon if needed
if ~isequal(polygon(1,:), polygon(end,:))
    polygon = [polygon; polygon(1,:)];
end

n_edges = size(polygon, 1) - 1;
if n_edges == 0
    minDist = 0;
    return;
end

% Vectorized distance calculation for all edges at once
p1 = polygon(1:n_edges, :);           % Start points
p2 = polygon(2:n_edges+1, :);         % End points

% Vectors from p1 to p2 for all edges
v = p2 - p1;                          % [n_edges x 2]

% Vectors from p1 to point for all edges
w = point - p1;                       % [n_edges x 2]

% Dot products
v_dot_v = sum(v .* v, 2);             % [n_edges x 1]
w_dot_v = sum(w .* v, 2);             % [n_edges x 1]

% Handle degenerate cases (zero-length edges)
valid_edges = v_dot_v > eps;
t = zeros(n_edges, 1);
t(valid_edges) = w_dot_v(valid_edges) ./ v_dot_v(valid_edges);

% Clamp t to [0,1]
t = max(0, min(1, t));

% Closest points on all segments
closest_points = p1 + t .* v;         % [n_edges x 2]

% Distances to all segments
distances = sqrt(sum((point - closest_points).^2, 2));

% Minimum distance
minDist = min(distances);

end