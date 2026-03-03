function test_progress_bar_object()
    n = 10000 ;
    pbo = progress_bar_object(n) ;
    for i = 1 : n ,
        pbo.update() ;
    end
end
