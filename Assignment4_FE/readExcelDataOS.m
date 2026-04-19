function [dates, rates] = readExcelDataOS(filename, formatData)
% readExcelData - Versione moderna ottimizzata per Mac M1/Apple Silicon
% Legge dati, calcola i tassi Mid (media tra Bid e Ask) e gestisce le date.

%% Dates from Excel (usiamo readcell per catturare testo e formati data)
% Settlement date
c_settlement = readcell(filename, 'Sheet', 1, 'Range', 'E8');
dates.settlementDate = parseDate(c_settlement{1}, formatData); 

% Dates relative to depos
c_depos = readcell(filename, 'Sheet', 1, 'Range', 'D11:D18');
dates.depos = parseDate(c_depos, formatData);

% Dates relative to futures: calc start & end (Q=Start, R=End)
c_futures = readcell(filename, 'Sheet', 1, 'Range', 'Q12:R20');
dates.futures = parseDate(c_futures, formatData);

% Date relative to swaps: expiry dates
c_swaps = readcell(filename, 'Sheet', 1, 'Range', 'D39:D88');
dates.swaps = parseDate(c_swaps, formatData);

%% Rates from Excel (usiamo readmatrix per i numeri)
% Calcoliamo in automatico il tasso MID (media tra Bid e Ask sulle righe)

% Depositi
tassi_depositi = readmatrix(filename, 'Sheet', 1, 'Range', 'E11:F18');
rates.depos = mean(tassi_depositi, 2) / 100;

% Futures
tassi_futures = readmatrix(filename, 'Sheet', 1, 'Range', 'E28:F36');
tassi_futures_mid = mean(tassi_futures, 2);
% Il file Excel dà il prezzo (es. 95). Il tasso implicito è 100 - Prezzo.
rates.futures = (100 - tassi_futures_mid) / 100;

% Swaps
tassi_swaps = readmatrix(filename, 'Sheet', 1, 'Range', 'E39:F88');
rates.swaps = mean(tassi_swaps, 2) / 100;

end

% =========================================================================
% HELPER FUNCTION: Risolve il problema delle date lette su Mac
% =========================================================================
function dn = parseDate(c, formatData)
    % Se c è una singola cella o una matrice di celle, converte in datenum
    if iscell(c)
        dn = zeros(size(c));
        for i = 1:numel(c)
            val = c{i};
            if isdatetime(val)
                dn(i) = datenum(val); % Mac spesso legge come datetime automatico
            elseif ischar(val) || isstring(val)
                dn(i) = datenum(val, formatData); % Se è testuale, usa il formato
            elseif isnumeric(val)
                dn(i) = val + 693960; % Se è il numero seriale di Excel (1900-based)
            end
        end
    else
        % Caso singolo
        if isdatetime(c)
            dn = datenum(c);
        elseif ischar(c) || isstring(c)
            dn = datenum(c, formatData);
        elseif isnumeric(c)
            dn = c + 693960; 
        end
    end
end