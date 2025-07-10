load('/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/FlyBubbleRGB/locomotionGtACR1_24_RGB_EXT_VGLUT-GAL4_RigA_20210305T083721/protocol.mat');
% old = load('/groups/branson/home/robiea/Projects_data/Labeler_APT/cx_GMR_SS00020_CsChr_RigB_20150908T133237/protocol.mat');

% load(fullfile(expdir,dataloc_params.ledprotocolfilestr),'protocol')

if isfield(protocol,'Rintensity')
    RGBprotocol = protocol;
    clear protocol;
    % test if RGBprotocol has only one active color
    countactiveLEDs = [double(any(RGBprotocol.Rintensity));double(any(RGBprotocol.Gintensity));double(any(RGBprotocol.Bintensity))];
    
    % check that there is 1 and only 1 color LED used in protocol
    if sum(countactiveLEDs) == 0 
        error('ChR = 1 for LED protcol with no active LEDs')
    elseif sum(countactiveLEDs) > 1
        error('More than one active LED color in protocol. Not currently supported')
    end
    
    % call function that transforms new protocol to old protocol
    [protocol,ledcolor] = convertRGBprotocol2protocolformat(RGBprotocol,countactiveLEDs)
    
    
    
end

    