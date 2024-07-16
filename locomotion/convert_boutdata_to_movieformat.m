function [movieformat_perfly_limbboutdata] = convert_boutdata_to_movieformat(trx,perfly_limbboutdata)

% input: 
% (output of perfly_limbboutdata or restructedSwingStance in
% trajectory format)

% perfly_limbboutdata struct with perlimb field for each fly:
%  fields:
%     start_indices_stance: first frame of groundcontact 
%     end_indices_stance: first frame of not groundcontact after paired start_indices_stance
%     start_indices_swing
%     end_indices_swing


% output same struct with values in movie reference frame, 1 = firstframe of
% movie for all flies
%     start_indices_stance: first frame of groundcontact 
%     end_indices_stance: first frame of not groundcontact after paired start_indices_stance
%     start_indices_swing
%     end_indices_swing
state = {'swing','stance'};

movieformat_perfly_limbboutdata = struct;
for fly = 1:numel(perfly_limbboutdata)
    for limb = 1:numel(perfly_limbboutdata(fly).perlimb)
        for is = 1:numel(state)
            movieformat_perfly_limbboutdata(fly).perlimb(limb).(state{is}).start_indices = perfly_limbboutdata(fly).perlimb(limb).(state{is}).start_indices-trx(fly).off;
            movieformat_perfly_limbboutdata(fly).perlimb(limb).(state{is}).end_indices = perfly_limbboutdata(fly).perlimb(limb).(state{is}).end_indices-trx(fly).off;

        end
% 
%         movieformat_perfly_limbboutdata(fly).perlimb(limb).start_indices_stance = perfly_limbboutdata(fly).perlimb(limb).start_indices_stance-trx(fly).off;
%         movieformat_perfly_limbboutdata(fly).perlimb(limb).end_indices_stance = perfly_limbboutdata(fly).perlimb(limb).end_indices_stance-trx(fly).off;
%         movieformat_perfly_limbboutdata(fly).perlimb(limb).start_indices_swing = perfly_limbboutdata(fly).perlimb(limb).start_indices_swing-trx(fly).off;
%         movieformat_perfly_limbboutdata(fly).perlimb(limb).end_indices_swing = perfly_limbboutdata(fly).perlimb(limb).end_indices_swing-trx(fly).off;
    end
end