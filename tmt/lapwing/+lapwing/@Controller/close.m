function close(self)

% fpj=get(handle(self.figure_h),'JavaFrame');
% jw=fpj.fHG1Client.getWindow;
% jcb=handle(jw,'CallbackProperties');
% set(jcb,'WindowGainedFocusCallback',[]);  % have to do this, or else
%                                           % the MVC objects hang around
% clear fpj jw jcb;

delete(self.figure_h) ;
self.figure_h = [] ;

end
