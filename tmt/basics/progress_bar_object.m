classdef progress_bar_object < handle
    properties
        n_
        i_
        percent_as_displayed_last_
        did_print_at_least_one_line_
        did_print_final_newline_
        data_queue_
    end
    
    methods
        function self = progress_bar_object(n)
            self.n_ = n ;
            self.i_ = 0 ;
            self.percent_as_displayed_last_ = [] ;
            self.did_print_at_least_one_line_ = false ;
            self.did_print_final_newline_ = false ;
            self.data_queue_ = parallel.pool.DataQueue ;
            self.data_queue_.afterEach(@(di)(self.update_helper_(di))) ;
        end
        
        function delete(self)
            if ~isempty(self.data_queue_) && isvalid(self.data_queue_) ,
                delete(self.data_queue_) ;
                self.data_queue_ = [] ;
            end
        end
        
        function update(self, varargin)
            % Should be called from within the parfor loop
            if isempty(varargin) ,
                di = 1 ;
            else
                di = varargin{1} ;
            end
            self.data_queue_.send(di) ;
        end
        
        function update_helper_(self, di)
            % Should not be called by clients
            if self.did_print_final_newline_ ,
                return
            end
            self.i_ = min(self.i_ + di, self.n_) ;
            i = self.i_ ;
            n = self.n_ ;
            percent = fif(n==0, 100, 100*(i/n)) ;
            percent_as_displayed = round(percent*10)/10 ;    
            if ~isequal(percent_as_displayed, self.percent_as_displayed_last_) ,
                if self.did_print_at_least_one_line_ ,
                    delete_bar = repmat('\b', [1 1+50+1+2+4+1]) ;
                    fprintf(delete_bar) ;
                end
                bar = repmat('*', [round(percent/2) 1]) ;
                fprintf('[%-50s]: %4.1f%%', bar, percent_as_displayed) ;
                self.did_print_at_least_one_line_ = true ;
            end
            if i==n ,
                if ~self.did_print_final_newline_ ,
                    fprintf('\n') ;
                    self.did_print_final_newline_ = true ;
                end
            end
            self.percent_as_displayed_last_ = percent_as_displayed ;
        end
    end
end
 