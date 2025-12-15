function pci = calculate_pci(U)
% Menghitung Partition Coefficient Index (PCI)
% Input: U (matriks partisi c x n)
% Output: pci (nilai PCI, range 0-1, semakin tinggi semakin baik)
% Formula: PCI = (1/n) * sum(sum(U^2))

    [c, n] = size(U);

    pci = sum(sum(U .^ 2)) / n;

    fprintf('\nPartition Coefficient Index (PCI)\n');
    fprintf('Formula: PCI = (1/n) * sum(sum(U^2))\n');
    fprintf('Jumlah cluster: %d, Jumlah data: %d\n', c, n);
    fprintf('Nilai PCI: %.6f\n', pci);
    fprintf('\nInterpretasi:\n');
    fprintf('  PCI = 1.0    : Partisi crisp sempurna\n');
    fprintf('  PCI = %.4f : Partisi fuzzy penuh (1/c)\n', 1/c);
    fprintf('  PCI tinggi = cluster lebih terdefinisi\n');

end
