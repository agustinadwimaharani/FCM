function si = calculate_silhouette_index(X, U)
    % X : data (n x d)
    % U : membership matrix FCM (c x n)

    % Tentukan cluster untuk setiap data (hard clustering)
    [~, labels] = max(U);

    n = size(X, 1);
    D = zeros(n);


    % 1. Hitung jarak 
    for i = 1:n
        for j = i+1:n
            dist = sqrt(sum((X(i,:) - X(j,:)).^2));
            D(i,j) = dist;
            D(j,i) = dist;
        end
    end

    % 2. Hitung nilai silhouette
    s = zeros(n,1);
    c = max(labels);

    for i = 1:n
        % Cluster data ke-i
        cluster_i = labels(i);

        % Hitung a(i): rata-rata jarak dengan cluster sendiri 
        same_cluster = find(labels == cluster_i);
        same_cluster(same_cluster == i) = [];  % kecuali dirinya sendiri

        if isempty(same_cluster)
            a = 0;
        else
            a = mean(D(i, same_cluster));
        end

        % Hitung b(i): jarak minimum ke cluster lain
        b = inf;

        for cluster_j = 1:c
            if cluster_j == cluster_i
                continue;
            end

            other_cluster = find(labels == cluster_j);

            if isempty(other_cluster)
                continue;
            end

            bj = mean(D(i, other_cluster));

            if bj < b
                b = bj;
            end
        end

        % Nilai silhouette s(i) 
        s(i) = (b - a) / max(a, b);
    end

  
    % 3. Silhouette Index keseluruhan
    si = mean(s);
end
