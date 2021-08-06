function flydisco_pipeline_error(stage, msgs)
    header = sprintf('Stage %s failed:', stage) ;
    indented_messages = cellfun(@(str)(sprintf('  %s', str)), msgs, 'UniformOutput', false) ;
    all_msgs = [ {header} , indented_messages ] ;
    final_error_message = sprintf('%s\n', all_msgs{:}) ;
    error('flydisco:pipeline_error', final_error_message) ;  %#ok<SPERR>
end
