function Granger_graph(sig,nodeNames,title_name)

sig_graph = sig;
sig_graph(logical(eye(size(sig_graph)))) = 0;

Graph = digraph(sig_graph',nodeNames);

% Plot the graph
h = plot(Graph, 'Layout', 'circle');

% Customize node appearance
h.MarkerSize = 15;             % Increase node size
h.NodeColor = 'k';             % Change node color to blue
h.ArrowSize = 15;              % Increase arrow size
h.EdgeColor = 'r';             % Edge color to red

% Get node positions
X = h.XData;
Y = h.YData;

% Remove default node labels to place custom ones
h.NodeLabel = '';

% Define offsets for label positioning
xOffset = 0.27;  % Horizontal offset
yOffset = 0.27;  % Vertical offset

% Add custom labels with offset positions
for k = 1:length(X)
    text(X(k) + xOffset, Y(k) + yOffset, nodeNames{k}, 'FontSize', 15, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'cap');
end

title(title_name);

end