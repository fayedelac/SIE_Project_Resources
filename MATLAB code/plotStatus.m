%% plotStatus.m file contains all necessary functions to plot the status
%% the first function calls all the other functions

function plotStatus(network_data, nodeId, statusObs, statusHierarchy, thresholdValue, thresholdUncertainty)
    % get Dates
    dates = getDates(network_data);

	ids = thresholdUncertainty < 20;
    ids = ones(3, 1, 'logical'); % AJOUT
	dates = dates(ids);
	statusObs = statusObs(ids, :);
	statusHierarchy = statusHierarchy(ids, :);
	thresholdValue = thresholdValue(ids, :);

	N = size(statusHierarchy, 2); % number of nodes
	T = size(statusHierarchy, 1); % number of timesteps

	% preparare threshold plotting
	xHie = [thresholdValue'; thresholdValue'] + 0.5;
	yHie = [0:T-1; 1:T] + 0.5;

	% prepare confusion matrix
    % 1: wet, -1: dry, 0: no observation
	C = NaN(T, N);
	C(statusObs == 1   & statusHierarchy == 1) = 1; % true positive
	C(statusObs == 0   & statusHierarchy == 1) = 2; % positive, not observed
	C(statusObs == -1  & statusHierarchy == 1) = 3; % false positive
	C(statusObs == 1   & statusHierarchy == -1) = 4; % false negative
	C(statusObs == 0   & statusHierarchy == -1) = 5; % negative, not observed
	C(statusObs == -1  & statusHierarchy == -1) = 6; % true negative

    % Print out the nodes that have false positive or false negative values
    % Find indices of false positives (3) and false negatives (4) in the confusion matrix
    falsePositive = C == 3;
    falseNegative = C == 4;
    falsePositiveColumns = any(falsePositive == 1, 1);
    falseNegativeColumns = any(falseNegative == 1, 1);
    falsePositiveIndices = find(falsePositiveColumns);
    falseNegativeIndices = find(falseNegativeColumns);
    falsePositiveNodeId = nodeId(falsePositiveIndices);
    falseNegativeNodeId = nodeId(falseNegativeIndices);
    %Check if there are any false positives
    if ~isempty(falsePositiveNodeId)
        fprintf('Nodes with False Positives:\n');
        disp(falsePositiveNodeId);
    else
        fprintf('No False Positives found.\n');
    end  
    % Check if there are any false negatives
    if ~isempty(falseNegativeNodeId)
        fprintf('Nodes with False Negatives:\n');
        disp(falseNegativeNodeId);
    else
        fprintf('No False Negatives found.\n');
    end

    % flip the matrix to display as in the paper of Durighetto et al. (2023)
    C = C.';
    C = flip(C);

	% figure
    f = figure;
    f.Position = [100 100 800 900];

	% plot timeseries
	imagesc(C);
    hold on;
    % plot line
    plot(xHie(:), yHie(:), '-', 'Color', 'k', 'LineWidth', 1.5);
    
    % x axis
	xticks(1:numel(dates));
	xticklabels(string(datestr(dates, 'dd.mm.yy')));
	xlabel('Observation dates', 'FontSize', 14);
	xlim([0.5, numel(dates)+0.5]);

    % y axis
    set(gca, 'YTickLabel', []); % We hide the labels 
	ylabel('Ranking in hierarchy graph', 'FontSize', 14, 'Position', [-0.2, numel(nodeId)/2, 0]);
	ylim([0.5, numel(nodeId)+0.5]);

	% colors and legend
    colormap([[0 0.5 1]; [0.6 0.8 1]; [1 0 0]; [0 0 1]; [1, 0.67, 0.5]; [0.91, 0.36, 0.09]]);
    cb = colorbar();
	cb.Ticks = 1+(5/6/2) :5/6: 6-(5/6/2);
	cb.TickLabels = ["TP", "P", "FP", "FN", "N", "TN"];
	cb.TickLength = 0;	

    % Define the arrow coordinates for a vertical arrow
    arrowX = [0.1 0.1];
    arrowY = [0.12 0.92]; % Switched the arrow direction
    
    % Create the arrow annotation
    arrow = annotation('arrow', arrowX, arrowY, 'LineWidth', 10, 'Color', [0.70, 0.90, 1.00]);
    arrow.HeadStyle = 'plain';
    arrow.HeadWidth = 19;
    arrow.HeadLength = 30;
    
    % Create an axes for the arrow
    ax1 = axes('Position', [0.01 0.1 0.2 0.8]); % Moved more to the left
    axis(ax1, 'off');
    
    % Add labels inside the arrow
    text(0.3, 0.99, 'Last', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', '#808080', 'Rotation', 90);
    text(0.3, 0.05, 'First', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', '#808080', 'Rotation', 90);
end

function dates = getDates(network_data)

    % retrieve the field names
    field_names = fieldnames(network_data);
    % the dates variables are found after the three first variables and before the
    % three last variables
    variable_dates = field_names(4:length(field_names)-3);
    % transform the date variable into a valid variable
    dates_string = cellfun(@(x) x(5:12), variable_dates, 'UniformOutput', false);

    dateFormat = 'dd_MM_yy';
    dates = datetime(dates_string, 'InputFormat', dateFormat);

end