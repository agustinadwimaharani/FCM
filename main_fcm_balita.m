% Clustering Status Gizi Balita - Fuzzy C-Means dengan PCI
% Clustering balita berdasarkan JK, Usia, Berat, Tinggi
% Evaluasi 3 skenario: 2, 3, dan 4 cluster

clear; clc; close all;
fprintf('\nCLUSTERING STATUS GIZI BALITA - FUZZY C-MEANS\n\n');

%% LOAD DATA
fprintf('Load Dataset\n');
data_path = 'DataBalita.xlsx';
fprintf('Membaca file: %s\n', data_path);

data = readtable(data_path,'Sheet',1);
data.Properties.VariableNames(1:6) = ...
    {'No','Nama','JK','Usia_Bulan','Berat','Tinggi'};
data = rmmissing(data);

fprintf('Data berhasil dimuat: %d balita\n\n', height(data));

%% PREPROCESSING 
[X_normalized, X_encoded, params] = preprocess_data(data);

%% FUZZY C-MEANS
cluster_scenarios = [2 3 4];
results = struct();

for idx = 1:length(cluster_scenarios)
    c = cluster_scenarios(idx);

    fprintf('\n--- SKENARIO %d CLUSTER ---\n', c);
    [U,V,J_history,iteration_data] = ...
        fcm_detailed(X_normalized,c,2,100,1e-5);

    results(idx).c = c;
    results(idx).U = U;
    results(idx).V = V;
    results(idx).J_history = J_history;
    results(idx).iteration_data = iteration_data;

    [~,lbl] = max(U,[],1);
    results(idx).cluster_labels = lbl';
    results(idx).pci = calculate_pci(U);
    results(idx).si  = calculate_silhouette_index(X_normalized,U);
end

%% EVALUASI PCI 
pci_values = [results.pci];
[best_pci,best_idx] = max(pci_values);

fprintf('\nEVALUASI PCI\n');
for i=1:length(results)
    fprintf('%d Cluster : %.6f\n',results(i).c,results(i).pci);
end
% fprintf('Cluster Optimal: %d (PCI %.6f)\n\n',results(best_idx).c,best_pci);

%% EVALUASI SI 
si_values = [results.si];

fprintf('EVALUASI SI\n');
for i=1:length(results)
    fprintf('%d Cluster : %.6f\n',results(i).c,results(i).si);
end


%% SIMPAN HASIL 
output_dir = 'output';
if ~exist(output_dir,'dir'), mkdir(output_dir); end

%% SIMPAN MATRIKS PARTISI 
m = 2;
for idx = 1:length(results)
    c = results(idx).c;

    U_last = results(idx).iteration_data{end}.U;
    U_last = U_last ./ sum(U_last,1);
    U_pow  = U_last.^m;

    T_U = array2table(U_last','VariableNames', ...
        arrayfun(@(x)sprintf('C%d',x),1:c,'UniformOutput',false));
    T_Um = array2table(U_pow','VariableNames', ...
        arrayfun(@(x)sprintf('C%d',x),1:c,'UniformOutput',false));

    writetable(T_U, fullfile(output_dir,...
        sprintf('Matriks_Partisi_%dCluster.xlsx',c)));
    writetable(T_Um, fullfile(output_dir,...
        sprintf('Matriks_Partisi_Pangkat_%dCluster.xlsx',c)));
end


%% SIMPAN FUNGSI OBJEKTIF (ITERASI 1 & ITERASI TERAKHIR PER SKENARIO)

for idx = 1:length(results)
    c = results(idx).c;

    iter1 = results(idx).iteration_data{1};
    iterLast = results(idx).iteration_data{end};

    % ===== ITERASI 1 =====
    J_detail_1 = iter1.J_detail;
    J_total_1  = iter1.J_total_per_data;
    n = size(J_detail_1,2);

    T_J1 = table((1:n)','VariableNames',{'Data'});
    for i = 1:c
        T_J1.(sprintf('C%d',i)) = J_detail_1(i,:)';
    end
    T_J1.Total = J_total_1';

    % ===== ITERASI TERAKHIR =====
    J_detail_L = iterLast.J_detail;
    J_total_L  = iterLast.J_total_per_data;

    T_JL = table((1:n)','VariableNames',{'Data'});
    for i = 1:c
        T_JL.(sprintf('C%d',i)) = J_detail_L(i,:)';
    end
    T_JL.Total = J_total_L';

    % ===== SIMPAN DALAM 1 FILE, 2 SHEET =====
    filename = fullfile(output_dir, ...
        sprintf('Fungsi_Objektif_%dCluster.xlsx', c));

    writetable(T_J1, filename, 'Sheet', 'Iterasi_1');
    writetable(T_JL, filename, 'Sheet', 'Iterasi_Terakhir');

    fprintf('Saved: Fungsi_Objektif_%dCluster.xlsx (Iterasi 1 & Terakhir)\n', c);
end

%% SIMPAN MATRIKS PERUBAHAN PER DATA (ITERASI 1 & ITERASI TERAKHIR)

for idx = 1:length(results)
    c = results(idx).c;

    iter1    = results(idx).iteration_data{1};
    iterLast = results(idx).iteration_data{end};

    % === Ambil Matriks Perubahan per Data ===
    DeltaU_1 = iter1.delta_U_matrix;    % c x n
    DeltaU_L = iterLast.delta_U_matrix; % c x n

    % === Transpose agar per baris = data (lebih rapi di Excel) ===
    DeltaU_1 = DeltaU_1';
    DeltaU_L = DeltaU_L';

    % === Konversi ke tabel ===
    T_DU1 = array2table(DeltaU_1, ...
        'VariableNames', arrayfun(@(x)sprintf('C%d',x),1:c,'UniformOutput',false));

    T_DUL = array2table(DeltaU_L, ...
        'VariableNames', arrayfun(@(x)sprintf('C%d',x),1:c,'UniformOutput',false));

    % === Simpan ke Excel (1 file, 2 sheet) ===
    filename = fullfile(output_dir, ...
        sprintf('Matriks_Perubahan_Per_Data_%dCluster.xlsx', c));

    writetable(T_DU1, filename, 'Sheet', 'Iterasi_1');
    writetable(T_DUL, filename, 'Sheet', 'Iterasi_Terakhir');

    fprintf('Saved: Matriks_Perubahan_Per_Data_%dCluster.xlsx\n', c);
end


%% Visualisasi
fprintf('Membuat Visualisasi\n');

for idx = 1:length(cluster_scenarios)
    V_denorm = zeros(size(results(idx).V));
    for i = 1:4
        V_denorm(:, i) = results(idx).V(:, i) * params.range(i) + params.min_vals(i);
    end
    results(idx).V_denorm = V_denorm;
end

for idx = 1:length(cluster_scenarios)
    c = results(idx).c;
    cluster_labels = results(idx).cluster_labels;
    V_denorm = results(idx).V_denorm;

    figure('Name', sprintf('%d Cluster - PCI %.4f', c, results(idx).pci), ...
           'Position', [100 + (idx-1)*50, 100 + (idx-1)*50, 1200, 400]);

    % Plot 1: Berat vs Tinggi
    subplot(1, 3, 1);
    colors = lines(c);
    for i = 1:c
        idx_cluster = find(cluster_labels == i);
        scatter(X_encoded(idx_cluster, 3), X_encoded(idx_cluster, 4), 50, ...
                colors(i, :), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        hold on;
    end
    scatter(V_denorm(:, 3), V_denorm(:, 4), 200, colors, 'p', ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    xlabel('Berat Badan (kg)', 'FontWeight', 'bold');
    ylabel('Tinggi Badan (cm)', 'FontWeight', 'bold');
    title(sprintf('Berat vs Tinggi (%d Cluster)', c), 'FontWeight', 'bold');
    legend_entries = cell(c*2, 1);
    for i = 1:c
        legend_entries{i} = sprintf('Cluster %d', i);
        legend_entries{c+i} = sprintf('Centroid %d', i);
    end
    legend(legend_entries, 'Location', 'best');
    grid on;
    hold off;

    % Plot 2: 3D
    subplot(1, 3, 2);
    for i = 1:c
        idx_cluster = find(cluster_labels == i);
        scatter3(X_encoded(idx_cluster, 2), X_encoded(idx_cluster, 3), ...
                 X_encoded(idx_cluster, 4), 50, colors(i, :), 'filled', ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        hold on;
    end
    scatter3(V_denorm(:, 2), V_denorm(:, 3), V_denorm(:, 4), 200, ...
             colors, 'p', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    xlabel('Usia (Bulan)', 'FontWeight', 'bold');
    ylabel('Berat (kg)', 'FontWeight', 'bold');
    zlabel('Tinggi (cm)', 'FontWeight', 'bold');
    title(sprintf('3D: Usia-Berat-Tinggi (%d Cluster)', c), 'FontWeight', 'bold');
    legend(legend_entries, 'Location', 'best');
    grid on;
    view(45, 30);
    hold off;

    % Plot 3: Konvergensi
    subplot(1, 3, 3);
    plot(1:length(results(idx).J_history), results(idx).J_history, ...
         'o-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    xlabel('Iterasi', 'FontWeight', 'bold');
    ylabel('Fungsi Objektif (J)', 'FontWeight', 'bold');
    title(sprintf('Konvergensi FCM (%d Cluster)', c), 'FontWeight', 'bold');
    grid on;

    sgtitle(sprintf('Hasil Clustering %d Cluster (PCI = %.4f)', c, results(idx).pci), ...
            'FontWeight', 'bold', 'FontSize', 14);
end

% Plot Perbandingan PCI
figure('Name', 'Perbandingan PCI', 'Position', [200, 200, 800, 500]);
bar(cluster_scenarios, pci_values, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on;
plot(cluster_scenarios, pci_values, 'ro-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Jumlah Cluster', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Partition Coefficient Index (PCI)', 'FontWeight', 'bold', 'FontSize', 12);
title('Perbandingan PCI untuk Berbagai Jumlah Cluster', 'FontWeight', 'bold', 'FontSize', 14);
xticks(cluster_scenarios);
grid on;
ylim([min(pci_values)*0.95, max(pci_values)*1.05]);
for i = 1:length(cluster_scenarios)
    text(cluster_scenarios(i), pci_values(i) + 0.005, sprintf('%.4f', pci_values(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
end
hold off;


%% Plot Perbandingan SI
figure('Name', 'Perbandingan Silhouette Index', ...
       'Position', [250, 250, 800, 500]);

bar(cluster_scenarios, si_values, ...
    'FaceColor', [0.4 0.8 0.4], ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.5);
hold on;

plot(cluster_scenarios, si_values, 'o-', ...
     'LineWidth', 2, ...
     'MarkerSize', 10, ...
     'MarkerFaceColor', 'g');

xlabel('Jumlah Cluster', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Silhouette Index (SI)', 'FontWeight', 'bold', 'FontSize', 12);
title('Perbandingan Silhouette Index untuk Berbagai Jumlah Cluster', ...
      'FontWeight', 'bold', 'FontSize', 14);

xticks(cluster_scenarios);
grid on;

ylim([min(si_values)*0.95, max(si_values)*1.05]);

for i = 1:length(cluster_scenarios)
    text(cluster_scenarios(i), si_values(i) + 0.005, ...
         sprintf('%.4f', si_values(i)), ...
         'HorizontalAlignment', 'center', ...
         'FontWeight', 'bold', ...
         'FontSize', 11);
end
hold off;

[best_si, best_si_idx] = max(si_values);

fprintf('\nCluster Optimal Berdasarkan SI : %d (SI %.6f)\n', ...
    results(best_si_idx).c, best_si);


%% Simpan Hasil
fprintf('Menyimpan Hasil ke Folder Output\n');

output_dir = 'output';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Simpan figure
fig_list = findall(0, 'Type', 'figure');
for i = 1:length(fig_list)
    fig_name = get(fig_list(i), 'Name');
    if isempty(fig_name)
        fig_name = sprintf('Figure_%d', i);
    end
    fig_name = strrep(fig_name, ' ', '_');
    fig_name = strrep(fig_name, ':', '');
    fig_name = strrep(fig_name, '-', '');

    saveas(fig_list(i), fullfile(output_dir, [fig_name, '.png']));
    fprintf('  Saved: %s.png\n', fig_name);
end

% Simpan hasil clustering
for idx = 1:length(cluster_scenarios)
    c = results(idx).c;

    output_table = table();
    output_table.No = (1:height(data))';
    output_table.Nama = data.Nama;
    output_table.JK = data.JK;
    output_table.Usia_Bulan = data.Usia_Bulan;
    output_table.Berat = data.Berat;
    output_table.Tinggi = data.Tinggi;

    for i = 1:c
        output_table.(sprintf('Membership_C%d', i)) = results(idx).U(i, :)';
    end

    output_table.Cluster = results(idx).cluster_labels;

    filename = sprintf('hasil_clustering_%d_cluster.xlsx', c);
    writetable(output_table, fullfile(output_dir, filename));
    fprintf('  Saved: %s\n', filename);
end

% Simpan PCI
pci_table = table(cluster_scenarios', pci_values', 'VariableNames', {'Jumlah_Cluster', 'PCI'});
writetable(pci_table, fullfile(output_dir, 'perbandingan_pci.xlsx'));
fprintf('  Saved: perbandingan_pci.xlsx\n');



% Simpan workspace
save(fullfile(output_dir, 'workspace_fcm.mat'));
fprintf('  Saved: workspace_fcm.mat\n');



fprintf('\nSELESAI! Semua hasil disimpan di folder: %s\n\n', output_dir);

%% RINGKASAN 
fprintf('\nRINGKASAN HASIL\n');
fprintf('Dataset : %d Balita\n',height(data));
fprintf('Metode  : Fuzzy C-Means\n');
fprintf('Cluster Optimal : %d (PCI %.6f)\n',...
    results(best_idx).c,best_pci);
