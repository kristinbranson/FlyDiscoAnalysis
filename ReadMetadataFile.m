function metadata = ReadMetadataFile(metadatafile)

skip_prefix_keywords = {'experiment','session','flies','apparatus','rearing','handling','environment'};

DOMnode= xmlread(metadatafile);
n = DOMnode.getDocumentElement();
metadata = struct;
parse(n,'');


  function parse(n,prefix)
    
    % add attributes
    
    name = char(n.getNodeName());
    a = n.getAttributes();
    if ~ismember(name,skip_prefix_keywords),
      if ~isempty(prefix),
        prefix = [prefix,'__'];
      end
      prefix = [prefix,char(n.getNodeName())];
      value = n.getNodeValue();
      if isempty(a) || ~isempty(value),
        metadata.(prefix) = value;
      end
    end
    
    for i = 0:a.getLength()-1,
      name = char(a.item(i).getName());
      if ~isempty(prefix),
        name = [prefix,'__',name]; %#ok<AGROW>
      end
      value = char(a.item(i).getValue());
      if isfield(metadata,name),
        warning('Overwriting file %s\n',name);
      end
      metadata.(name) = value;
    end
    
    cs = n.getChildNodes();
    for i = 0:cs.getLength()-1,
      c = cs.item(i);
      if c.getNodeType() == 1,
        parse(c,prefix);
      end
    end
    
  end

end