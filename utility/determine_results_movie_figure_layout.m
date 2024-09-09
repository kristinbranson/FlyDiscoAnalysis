function [final_nzoomr, final_nzoomc, final_figpos] = ...
        determine_results_movie_figure_layout(ctraxresultsmovie_params, trx)

    if ischar(ctraxresultsmovie_params.nzoomr) || ischar(ctraxresultsmovie_params.nzoomc),
        
        % figure out number of flies
        %load(trxfile,'trx');
        
        firstframe = min([trx.firstframe]);
        endframe = max([trx.endframe]);
        trxnframes = endframe-firstframe+1;
        nflies = zeros(1,trxnframes);
        for i = 1:numel(trx),
            j0 = trx(i).firstframe-firstframe+1;
            j1 = trx(i).endframe-firstframe+1;
            nflies(j0:j1) = nflies(j0:j1)+1;
        end
        mediannflies = median(nflies);
        
        if isnumeric(ctraxresultsmovie_params.nzoomr),
            nzoomr = ctraxresultsmovie_params.nzoomr;
            nzoomc = round(mediannflies/nzoomr);
        elseif isnumeric(ctraxresultsmovie_params.nzoomc),
            nzoomc = ctraxresultsmovie_params.nzoomc;
            nzoomr = round(mediannflies/nzoomc);
        else
            nzoomr = ceil(sqrt(mediannflies));
            nzoomc = round(mediannflies/nzoomr);
        end
        final_nzoomr = nzoomr;
        final_nzoomc = nzoomc;
        
        if iscell(ctraxresultsmovie_params.figpos),
            [readframe,~,moviefile_fid] = get_readframe_fcn(moviefile);
            if moviefile_fid < 0 ,
                error('Unable to read file %s', moviefile) ;
            end
            moviefile_cleaner = onCleanup(@()(fclose(moviefile_fid))) ;
            im = readframe(1);
            [nr,nc,~] = size(im);
            
            rowszoom = floor(nr/nzoomr);
            imsize = [nr,nc+rowszoom*nzoomc];
            figpos = str2double(ctraxresultsmovie_params.figpos);
            if isnan(figpos(3)),
                figpos(3) = figpos(4)*imsize(2)/imsize(1);
            elseif isnan(figpos(4)),
                figpos(4) = figpos(3)*imsize(1)/imsize(2);
            end
            final_figpos = figpos;
        else
            final_figpos = ctraxresultsmovie_params.figpos ;
        end        
    else
        final_nzoomr = ctraxresultsmovie_params.nzoomr ;
        final_nzoomc = ctraxresultsmovie_params.nzoomc ;
        final_figpos = ctraxresultsmovie_params.figpos ;
    end
    
end