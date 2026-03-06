function m = DecodeStatName(s)

m = regexp(s,'^(?<field>.*)_fly(?<flycondition>.*)_frame(?<framecondition>.*)$','names','once');