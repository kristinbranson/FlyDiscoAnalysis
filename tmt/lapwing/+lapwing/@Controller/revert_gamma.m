function revert_gamma(self)
    self.model.revert_gamma() ;
    set(self.figure_h,'colormap',self.model.cmap);
end
