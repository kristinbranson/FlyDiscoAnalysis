function [edit_pad_width, edit_pad_height, edit_horizontal_alignment, edit_bg_color] = get_platform_dependent_quantities()
    if ismac ,
        edit_pad_width=14;
        edit_pad_height=7;
        % setting horizontal alignment to 'right' on mac looks ugly
        % seems like a matlab bug
        edit_horizontal_alignment='center';
        % can't set this to white on mac -- background overflows the box
        % again, seems like a matlab bug.
        % setting it to something else looks even worse.
        edit_bg_color=[1 1 1];
    else
        edit_pad_width=10;
        edit_pad_height=2;
        edit_horizontal_alignment='right';
        edit_bg_color=[1 1 1];
    end
end
