%AR 20240602 
function [perlimbboutdata] = limbSwingStanceStep(groundcontactdata)
%input ground contact in trajectory data format (first frame of trajectory = 1)
%
% input groundcontact data flies x 6 leg tips
% (from params files [12,13,14,15,16,17] )
% output perfly_limbboutdata struct with perlimb field for each fly:
%  fields:
%     start_indices_stance: first frame of groundcontact 
%     end_indices_stance: first frame of not groundcontact after paired start_indices_stance
%     start_indices_swing
%     end_indices_swing
perlimbboutdata = struct;
for ifly = 1:numel(groundcontactdata)
    % 6 legs x frames
    currdata = groundcontactdata{ifly};

    perlimb = struct;
    for ilegs = 1:size(currdata,1)
        currlimbdata = currdata(ilegs,:);
        % compute stance data - groundcontact stance = 1
        [start_indices_stance,end_indices_stance] = detect_bouts(currlimbdata);
        perlimb(ilegs).stance.start_indices = start_indices_stance(1:end-1);
        perlimb(ilegs).stance.end_indices = end_indices_stance(1:end-1);

        % compute swing data - groundcontact swing = 0
%         [start_indices_swing,end_indices_swing] = detect_bouts(~currlimbdata);
%         perlimb(ilegs).swing.start_indices = start_indices_swing;
%         perlimb(ilegs).swing.end_indices = end_indices_swing;

        % compute swing data from stance indices
        perlimb(ilegs).swing.start_indices = end_indices_stance(1:end-1);
        perlimb(ilegs).swing.end_indices = start_indices_stance(2:end);

        %compute step starts and ends
%         start_indices_step = start_indices_stance(1:end-1);
%         end_indices_step = start_indices_stance(2:end);
        perlimb(ilegs).step.start_indices = start_indices_stance(1:end-1);
        perlimb(ilegs).step.end_indices = start_indices_stance(2:end);
    end
    perlimbboutdata(ifly).perlimb = perlimb;
end