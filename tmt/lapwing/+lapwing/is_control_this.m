function ret=is_control_this(char,event)

% Tests whether the user pressed Control/Command-<this>, in a 
% cross-platform way.  char should be a 1x1 char array holding the
% _lowercase_ char in question.  Returns false if event holds a
% control-shift-<this> or anything like that: it has to be
% control/command-c, with no other modifiers.

ret=false;
if strcmp(event.Character,char) && length(event.Modifier)==1
  if ismac()
    if strcmp(event.Modifier,'command')
      ret=true;
    end
  else
    if strcmp(event.Modifier,'control')
      ret=true;
    end
  end
end

end
