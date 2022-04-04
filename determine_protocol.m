function protocol = determine_protocol(metadata, ledprotocol_file_contents)
    % Examine the metadata and the contents of the LED protocol file to determine
    % the protocol.
    if strcmp(metadata.assay,'FlyBubbleRGB') || strcmp(metadata.assay,'FlyBowlRGB') ,
        maybeRGBprotocol = ledprotocol_file_contents.protocol ;
        if isfield(maybeRGBprotocol,'Rintensity')
            RGBprotocol =  maybeRGBprotocol ;
            % test if RGBprotocol has only one active color
            is_active_from_LED_index = [ any(RGBprotocol.Rintensity) ; any(RGBprotocol.Gintensity) ; any(RGBprotocol.Bintensity) ] ;
            % check that there is 1 and only 1 color LED used in protocol
            active_LED_count = sum(double(is_active_from_LED_index)) ;
            if active_LED_count == 0 ,
                error('ChR = 1 for LED protcol with no active LEDs')
            elseif active_LED_count == 1 ,
                % do nothing, this is what we want
            else
                error('More than one active LED color in protocol. Not currently supported')
            end
            % call function that transforms new protocol to old protocol
            protocol = ConvertRGBprotocol2protocolformat(RGBprotocol, is_active_from_LED_index) ;
        else
            protocol = maybeRGBprotocol ;                    
        end
    else
        protocol = ledprotocol_file_contents.protocol ;
    end    
end
