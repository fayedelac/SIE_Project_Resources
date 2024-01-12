%% getHierarchy.m file contains all the necessary functions to perfrom the hierarchy reconstruction 
%% the first function calls all the other functions

function [network_obs_p, statusObs_p, observationMatrix_p, hierarchyMatrix_p, hierarchyGraph_p, order_p, nodeId_p, network_data_p, hierarchyAccuracy_p, statusHierarchy_p, thresholdValue_p, thresholdUncertainty_p, nodewiseAccuracy_p, timewiseAccuracy_p, totalAccuracy_p] = getHierarchy_p(network_data, network_obs, nodeId, percentage_values, len_p, k)
    % This function returns the hierarchy matrix 


    %% Initialize matrices containing variables for each percentage values
    % general variables
    network_obs_p = ones(size(network_obs, 1), size(network_obs, 2), len_p);
    statusObs_p = ones(size(network_obs, 2), size(network_obs, 1), len_p);
    observationMatrix_p = ones(size(network_obs, 1), size(network_obs, 1), len_p);
    hierarchyMatrix_p = ones(size(network_obs, 1), size(network_obs, 1), len_p);
    hierarchyGraph_p = cell(1, len_p);
    order_p = ones(1, size(network_obs, 1), len_p);
    nodeId_p = ones(size(network_obs, 1), 1, len_p);
    network_data_p = cell(1, len_p);
    % accuracy variables
    hierarchyAccuracy_p = zeros(1,1,len_p);
    statusHierarchy_p = zeros(size(network_obs, 2),size(network_obs, 1), len_p);
    thresholdValue_p = zeros(size(network_obs, 2), 1, len_p);
    thresholdUncertainty_p = zeros(size(network_obs, 2), 1, len_p);
    nodewiseAccuracy_p = zeros(1, size(network_obs, 1), len_p);
    timewiseAccuracy_p = zeros(size(network_obs, 2), 1, len_p);
    totalAccuracy_p = zeros(1,1,len_p);
    
    %% Temporary matrices, for every iteration of k
    % general variables
    network_obs_temp = ones(size(network_obs, 1), size(network_obs, 2), k);
    statusObs_temp = ones(size(network_obs, 2), size(network_obs, 1), k);
    observationMatrix_temp = ones(size(network_obs, 1), size(network_obs, 1), k);
    hierarchyMatrix_temp = ones(size(network_obs, 1), size(network_obs, 1), k);
    hierarchyGraph_temp = cell(1, k);
    order_temp = ones(1, size(network_obs, 1), k);
    nodeId_temp = ones(size(network_obs, 1), 1, k);
    network_data_temp = cell(1, k);
    % accuracy variables
    hierarchyAccuracy_temp = zeros(1,1,k);
    statusHierarchy_temp = zeros(size(network_obs, 2),size(network_obs, 1), k);
    thresholdValue_temp = zeros(size(network_obs, 2), 1, k);
    thresholdUncertainty_temp = zeros(size(network_obs, 2), 1, k);
    nodewiseAccuracy_temp = zeros(1, size(network_obs, 1), k);
    timewiseAccuracy_temp = zeros(size(network_obs, 2), 1, k);
    totalAccuracy_temp = zeros(1,1,k);


    % Loop for every different percentage values
    for i = 1:len_p
        fprintf('\nTotal accuracies for %0.f%% of available observations:\n', percentage_values(i)*100)

        % Cross-Validation: Loop for every k
        for j = 1:k

            % Get network observation matrix 
            network_obs_k = getNetworkObs(network_obs, percentage_values(i));

            % Create a matrix representing the observed states of each node on each date.
            % Each row corresponds to a day, and each column corresponds to a node.
            statusObs = getStatusMatrix(network_obs_k);

            % Create the observation graph 'O'. In this graph, an edge from node A to B exists
            % if there is an observation of wet A and dry B, suggesting A is more persistent than B.
            % The edge's weight is the number of such observations.
            observationMatrix = getObservationMatrix(statusObs);

            % Create the hierarchy graph 'H' by removing loops from the observation graph.
            % This graph represents the inferred hierarchy within the network.
            % breaks all cycles
            hierarchyMatrix = getHierarchyMatrix(observationMatrix);
            hierarchyMatrix = removeLoops(hierarchyMatrix);
            hierarchyGraph = digraph(hierarchyMatrix, network_data);

            %% Reorder graph topologically
            % toposort(G) returns the topological order of the nodes in G such that i < j for every edge (n(i),n(j)) in G. The directed graph G cannot have any cycles.
            order = toposort(hierarchyGraph);

            % Reorder variables
            network_data_k = network_data(order, :);
            nodeId_k = nodeId(order);
            statusObs = statusObs(:, order);
            observationMatrix = observationMatrix(order, order);    
            hierarchyMatrix = hierarchyMatrix(order, order);

            % Check hierarchy accuracy

            % Option 1: calculate the fraction of pairwise observations in accordance with the hierarchy. Since we reordered the matrix, these are all in the lower triangle.
            hierarchyAccuracy = sum(tril(observationMatrix), 'all') / sum(observationMatrix, 'all');

            % Option 2: use the hierarchy to reconstruct the status of the nodes, and calculate the corresponding accuracy in the usual way.
            [statusHierarchy, thresholdValue, thresholdUncertainty] = getStatusHierarchy(network_data_k, statusObs);

            nodewiseAccuracy = sum(statusHierarchy == statusObs) ./ sum(statusObs ~= 0);
            timewiseAccuracy = (statusHierarchy == statusObs) * network_data_k.LENGTH ./ ((statusObs ~= 0) * network_data_k.LENGTH);
            totalAccuracy = sum((statusHierarchy == statusObs) * network_data_k.LENGTH) / sum((statusObs ~= 0) * network_data_k.LENGTH);

            fprintf('Iteration #%0.f: %0.5f\n', j, totalAccuracy)
            
            % general variables
            network_obs_temp(:,:,j) = network_obs_k;
            statusObs_temp(:,:,j) = statusObs;
            observationMatrix_temp(:,:,j) = observationMatrix;
            hierarchyMatrix_temp(:,:,j) = hierarchyMatrix; 
            hierarchyGraph_temp{j} = hierarchyGraph;
            order_temp(:,:,j) = order;
            nodeId_temp(:,:,j) = nodeId_k;
            network_data_temp{j} = network_data_k;
            % accuracy variables
            hierarchyAccuracy_temp(:,:,j) = hierarchyAccuracy;
            statusHierarchy_temp(:,:,j) = statusHierarchy;
            thresholdValue_temp(:,:,j) = thresholdValue;
            thresholdUncertainty_temp(:,:,j) = thresholdUncertainty;
            nodewiseAccuracy_temp(:,:,j) = nodewiseAccuracy;
            timewiseAccuracy_temp(:,:,j) = timewiseAccuracy;
            totalAccuracy_temp(:,:,j) = totalAccuracy;

        end

    % Find iteration for which total accuracy closest to median total accuracy
    % Compute median accuracy
    median_accuracy = median(totalAccuracy_temp);
    fprintf('Median accuracy: %0.5f\n', median_accuracy)
    % Find the index of the value closest to median_accuracy
    [~, index] = min(abs(totalAccuracy_temp - median_accuracy));
    fprintf('Iteration #%0.f has the accuracy value the closest to median_accuracy: %0.5f\n', index, totalAccuracy_temp(:,:,index))
    % Save all necessary variables         
    network_obs_p(:,:, i) = network_obs_temp(:,:,index);
    observationMatrix_p(:,:, i) = observationMatrix_temp(:,:,index);
    hierarchyMatrix_p(:,:, i) = hierarchyMatrix_temp(:,:,index); 
    hierarchyGraph_p{i} = hierarchyGraph_temp{index};
    order_p(:,:, i) = order_temp(:,:,index);
    nodeId_p(:,:, i) = nodeId_temp(:,:,index);
    network_data_p{i} = network_data_temp{index};
    statusObs_p(:,:, i) = statusObs_temp(:,:,index);
    statusHierarchy_p(:,:, i) = statusHierarchy_temp(:,:,index);
    thresholdValue_p(:,:, i) = thresholdValue_temp(:,:,index);
    thresholdUncertainty_p(:,:, i) = thresholdUncertainty_temp(:,:,index);
    hierarchyAccuracy_p(:,:, i) = hierarchyAccuracy_temp(:,:,index);
    nodewiseAccuracy_p(:,:, i) = nodewiseAccuracy_temp(:,:,index);
    timewiseAccuracy_p(:,:, i) = timewiseAccuracy_temp(:,:,index);
    totalAccuracy_p(1,1,i) = totalAccuracy_temp(:,:,index);

    end

end

% ADDITIONAL FUNCTION - EXPERIMENTAL
% THIS FUNCTION IS USED ONLY FOR ANALYSIS, NOT FOR REAL STUDY
% network_obs returns the  observation arrays for percentage of
% observed nodes in the network
function network_obs_p = getNetworkObs(network_obs, percentage_value)

    % Determine the total number of elements in our observation array
    num_elements = numel(network_obs);

    temp = network_obs;

    % Determine the percentage of values we want to set to zero
    p = 1 - percentage_value;

    % Determine the number of elements to set to 0
    num_zeros = round((p/1) * num_elements);

    % Generate random indices to set to 0
    indices_to_zero = randperm(num_elements, num_zeros);

    % Set the selected elements to 0
    temp(indices_to_zero) = 0;

    network_obs_p = temp;      

end

% statusObs gets data from the provided table and transforms it into matrix form.
% -1: dry state, 1: wet state, 0: no observation
function statusObs = getStatusMatrix(network_obs)
    % One row per day and one column per node.
    statusObs = network_obs.';	
end

% O(i,j) = #obs where i was dry (-1) and j was wet (1)
function O = getObservationMatrix(statusObs)
	N = size(statusObs, 2); % number of nodes

	O = NaN(N, N);
	for n1 = 1:N
		for n2 = 1:N
			O(n1, n2) = sum(statusObs(:, n1) == -1 & statusObs(:, n2) == 1);
		end
	end
end

% The hierarchy matrix represents the direct acyclic graph (DAG) by
% breaking all cycles(i,j) and keeping i where i is more persistent than j
function hierarchyMatrix = getHierarchyMatrix(observationMatrix)
	% net observations
	hierarchyMatrix = observationMatrix' - observationMatrix;
	hierarchyMatrix(hierarchyMatrix < 0) = 0;
end

function H = removeLoops(H)
	maxNumCycles = 10^2;
	siz = size(H);

	% remove self loops -> already done in getHierarchy matrix?
	H = H - H';
	H(H < 0) = 0;

	% remove cycles
	for maxCycleLength = 2:2:siz(1)
		%disp(strcat("Removing cycles with size ", string(maxCycleLength)));
		% dot = 0;

		while true
			Hgraph = digraph(H);
			loops = allcycles(Hgraph, 'MaxNumCycles', maxNumCycles, 'MaxCycleLength', maxCycleLength);

			% disp(loops);
			% pause();

			if isempty(loops)
				break;
			end

			H = rmLoops(H, loops, siz);

			% dot = dot + 1;
			% if dot == 100; dot = 1; fprintf('\n'); end
			% fprintf('.');
		end
	end


	function H = rmLoops(H, loops, siz)
		% for each loop
		for i = 1:numel(loops)
			% get link strengths
			loop = [loops{i} loops{i}(1)];

			rr = loop(1:end-1);
			cc = loop(2:end);
			ii = (cc-1)*siz(1) + rr;

			strength = H(ii);

			% remove weaker link
			% disp(' ');
			% disp(loop);
			% disp(strength);

			[~, weakerPos] = min(strength);
			H(ii(weakerPos)) = 0;
		end

		% pause();
	end
end

function [statusHierarchy, thresholds, uncertainty] = getStatusHierarchy(network_data, statusObs)
% getStatusHierarchy determines the network status that minimizes the error. 
% To assess this, the function computes the status for various thresholds. 
% When testing the threshold for a specific node, nodes preceding it in the hierarchy are categorized as "wet," while nodes succeeding it are classified as "dry." 
% This process is repeated for different thresholds to find the optimal network status.

  
	T = size(statusObs, 2); % number of nodes, i.e. possible thresholds 
	N = size(statusObs, 1); % number of surveys

	stretchLength = network_data.LENGTH; % contains each length of all nodes

	errorHie = NaN(N, T); % each row contains the error relative to survey N, for each possible threshold

    for threshold = 0:T % testing all possible thresholds
        statusTemp = zeros(N, T);
        statusTemp(:, 1:threshold) = 1; % Nodes preceding the node tested (threshold) are categorized as "wet" = 1
        statusTemp(:, threshold+1:end) = -1;  % Nodes succeeding the node tested (threshold) are categorized as "dry" = -1
        
        errorHie(:, threshold+1) = (statusTemp ~= statusObs) * stretchLength; %~= not equal to, dim = NxT * Tx1 = Nx1
        % The error is proportional to the length of the nodes
    end

	% find treshold range associated to minimum error for each survey n
	thresholdMin = NaN(N, 1);
	thresholdMax = NaN(N, 1);
    for n = 1:N
        % contains 1 if the error is equal to the minimum error
		isMinError = errorHie(n, :) == min(errorHie(n, :));

        % find(X) returns a vector containing the linear indices of each nonzero element in array X.
        % We test first and last values because there might be two threshold that have the minimum error
		thresholdMin(n) = find(isMinError, 1, "first") - 1; % direction = 'first', finds the first index corresponding to nonzero elements.
		thresholdMax(n) = find(isMinError, 1, "last") - 1; % direction = 'last', finds the last index corresponding to nonzero elements in X
    end

	thresholds = round((thresholdMin + thresholdMax) / 2);
	uncertainty = thresholdMax - thresholdMin;

	% apply selected threshold
	statusHierarchy = NaN(N, T);

	for n = 1:N
        % 1: wet state, -1: dry state
		statusHierarchy(n, :) = [ones(1, thresholds(n)), -ones(1, T-thresholds(n))];
	end
end
