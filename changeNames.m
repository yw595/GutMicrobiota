function changedNames = changeNames(names)

changedNames = {};
for k=1:length(names)
    if ~isempty(regexp(names{k},'Clostridium'))
        changedNames{end+1} = 'Clostridium';
    elseif ~isempty(regexp(names{k},'Escherichia/Shigella'))
        changedNames{end+1} = 'Escherichia';
        changedNames{end+1} = 'Shigella';
    elseif ~isempty(regexp(names{k},'TM7'))
        changedNames{end+1} = 'candidate division TM7';
    elseif ~isempty(regexp(names{k},'Lachnospiracea'))
        changedNames{end+1} = 'Lachnospiraceae';
    else
        changedNames{end+1} = names{k};
    end
    changedNames{end} = strrep(changedNames{end},' ','_');
end

end