function handle_focus_gained(self)

% update the command-key state
%fprintf('focus gained\n');
self.shift_depressed=false;  % this may or may not really be correct,
                             % but it's best to assume this, so that
                             % user doesn't automatically get a
                             % circle/square when normal-dragging in
                             % certain circumstances.
                             
end
