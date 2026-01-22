function [U, V, J_history, iteration_data] = fcm_detailed(X, c, m, max_iter, epsilon)
% Fuzzy C-Means dengan output detail setiap iterasi
% Command Window: 10 data AWAL
% Penyimpanan: SEMUA data disimpan di iteration_data

    if nargin < 3, m = 2; end
    if nargin < 4, max_iter = 100; end
    if nargin < 5, epsilon = 1e-5; end

    [n, d] = size(X);

    fprintf('\nFuzzy C-Means - %d Cluster\n', c);
    fprintf('Parameter: n=%d, d=%d, m=%.1f, max_iter=%d, epsilon=%.0e\n', ...
            n, d, m, max_iter, epsilon);

    % Inisialisasi U
    U = rand(c, n);
    U = U ./ sum(U, 1);

    J_history = zeros(max_iter, 1);
    iteration_data = cell(max_iter, 1);

    for iter = 1:max_iter
        fprintf('\n--- Iterasi %d ---\n', iter);

        %% Hitung U^m
        Um = U .^ m;

        %% Hitung centroid
        V = (Um * X) ./ sum(Um, 2);

        fprintf('\nCentroid:\n');
        fprintf('%-10s %-12s %-12s %-12s %-12s\n', ...
                'Cluster', 'JK', 'Usia', 'Berat', 'Tinggi');
        for i = 1:c
            fprintf('%-10d %-12.6f %-12.6f %-12.6f %-12.6f\n', ...
                    i, V(i,1), V(i,2), V(i,3), V(i,4));
        end

        %% FUNGSI OBJEKTIF 
        dist2 = zeros(c, n);
        J_detail = zeros(c, n);
        J_total_per_data = zeros(1, n);

        for i = 1:c
            for j = 1:n
                dist2(i,j) = norm(X(j,:) - V(i,:))^2;
                J_detail(i,j) = Um(i,j) * dist2(i,j);
            end
        end

        % Total fungsi objektif
        J = sum(J_detail, 'all');
        J_history(iter) = J;

        fprintf('\nFungsi Objektif (J): %.6f\n', J);

        for j = 1:n
            J_total_per_data(j) = sum(J_detail(:,j));
        end

        %% TAMPILKAN 10 DATA AWAL
        fprintf('\nFungsi Objektif Detail per Data (10 data awal):\n');
        fprintf('%-8s', 'Data');
        for i = 1:c
            fprintf('%-12s', sprintf('C%d', i));
        end
        fprintf('%-12s\n', 'Total');

        for j = 1:min(10,n)
            fprintf('%-8d', j);
            for i = 1:c
                fprintf('%-12.6f', J_detail(i,j));
            end
            fprintf('%-12.6f\n', J_total_per_data(j));
        end

        %%  UPDATE U 
        U_old = U;
        D = zeros(c, n);

        for i = 1:c
            for j = 1:n
                D(i,j) = norm(X(j,:) - V(i,:));
            end
        end
        D(D == 0) = eps;

        for i = 1:c
            for j = 1:n
                sum_term = 0;
                for k = 1:c
                    sum_term = sum_term + (D(i,j)/D(k,j))^(2/(m-1));
                end
                U(i,j) = 1 / sum_term;
            end
        end

        %% DELTA U 
        delta_U = norm(U - U_old);
        fprintf('Delta U (norm total): %.6e\n', delta_U);

        delta_U_matrix = abs(U - U_old);

        fprintf('\nMatriks Perubahan per Data (10 data awal):\n');
        fprintf('%-8s', 'Data');
        for i = 1:c
            fprintf('%-12s', sprintf('C%d', i));
        end
        fprintf('\n');

        for j = 1:min(10,n)
            fprintf('%-8d', j);
            for i = 1:c
                fprintf('%-12.6f', delta_U_matrix(i,j));
            end
            fprintf('\n');
        end

        %% SIMPAN DETAIL ITERASI 
        iteration_data{iter}.J_detail = J_detail;                 % SEMUA DATA
        iteration_data{iter}.J_total_per_data = J_total_per_data; % SEMUA DATA
        iteration_data{iter}.V = V;
        iteration_data{iter}.U = U;
        iteration_data{iter}.J = J;
        iteration_data{iter}.delta_U_matrix = delta_U_matrix;
        iteration_data{iter}.delta_U = delta_U;

        %% CEK KONVERGENSI 
        if delta_U < epsilon
            fprintf('\nKonvergensi pada iterasi %d\n', iter);
            J_history = J_history(1:iter);
            iteration_data = iteration_data(1:iter);
            break;
        end
    end

    %% HASIL AKHIR
    fprintf('\nHasil Akhir:\n');
    fprintf('Fungsi Objektif final: %.6f\n', J_history(end));

    [~, cluster_labels] = max(U, [], 1);
    fprintf('\nDistribusi Cluster:\n');
    for i = 1:c
        fprintf('Cluster %d: %d data\n', i, sum(cluster_labels == i));
    end
end
