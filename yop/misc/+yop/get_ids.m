function id = get_ids(node_cell)
id = cellfun(@(e) e.id, node_cell);
end