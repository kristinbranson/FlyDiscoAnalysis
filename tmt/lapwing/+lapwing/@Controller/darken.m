function darken(self)
    self.model.darken() ;
    set(self.figure_h,'Colormap',self.model.cmap) ;
end
