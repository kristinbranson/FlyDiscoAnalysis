function flydisco_pipeline_error(stage, msgs)
    all_msgs = [ {sprintf('Stage %s failed:', stage)} , msgs ] ;
    error('%s\n', all_msgs{:});    
end
