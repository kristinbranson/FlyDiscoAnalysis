function brighten(self)
    self.model.brighten() ;
    set(self.figure_h,'Colormap',self.model.cmap) ;
end
