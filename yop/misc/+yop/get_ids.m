function id = get_ids(node_cell)
id = cellfun(@(e) get_id(e), node_cell);
end