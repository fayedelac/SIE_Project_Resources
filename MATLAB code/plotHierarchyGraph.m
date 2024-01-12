%% plotHierarchyGraph contains the unique function necessary to plot the hierarchy graph

function plotHierarchyGraph(hierarchyGraph)
	HreducedGraph = transreduction(hierarchyGraph);

	figure('Position', [100, 100, 400, 400]);
	%plot(HreducedGraph, 'Layout', 'layered', 'Direction', 'right', 'NodeColor', [0 0.4470 0.7410], 'EdgeColor', [0.8, 0.8, 0.8], 'MarkerSize', 6);
    plot(HreducedGraph, 'Layout', 'layered', 'Direction', 'down', 'NodeColor', [0 0.4470 0.7410], 'EdgeColor', [0.8, 0.8, 0.8], 'MarkerSize', 6);

    % y axis
    set(gca, 'YTickLabel', []); % Hide the labels 

    % Define the arrow coordinates for a vertical arrow
    arrowX = [0.1 0.1];
    arrowY = [0.12 0.92]; 
    
    % Create the arrow annotation
    arrow = annotation('arrow', arrowX, arrowY, 'LineWidth', 7, 'Color', [0.70, 0.90, 1.00]);
    arrow.HeadStyle = 'plain';
    arrow.HeadWidth = 15;
    arrow.HeadLength = 26;
    
    % Create an axes for the arrow
    ax1 = axes('Position', [0.01 0.1 0.2 0.8]); 
    axis(ax1, 'off');
    
    % Add labels inside the arrow
    text(0.22, 0.95, 'Last', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', '#808080', 'Rotation', 90);
    text(0.22, 0.08, 'First', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', '#808080', 'Rotation', 90);
end