function id = get_ids(node_cell)
id = cellfun(@(e) ID(e), node_cell);
end