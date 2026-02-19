function dataTable = getPigeon_biasCorrections(dataTable, biasCorrectionFile)
%function dataTable = getPigeon_biasCorrections(dataTable, biasCorrectionFile)

load(biasCorrectionFile)

SNRs = abs(dataTable.snr);
Lsnr = [SNRs<0.4 SNRs>0.4&SNRs<0.8 SNRs>0.8]; % yeah, ick. sue me.
numSNRs = size(Lsnr,2);
numDTs = size(boundBiasCorrection, 1);
for dd = 1:numDTs
    if dd < numDTs
        Ld = dataTable.DT == dd;
    else
        Ld = dataTable.DT >= dd;
    end
    for mm = 1:numSNRs
        Ls = Ld & Lsnr(:,mm);
        if any(Ls)
            newBounds = max(0.005, ...
                boundBiasCorrection(dd,1,mm) + ...
                boundBiasCorrection(dd,2,mm) .* abs(dataTable.bound(Ls)));

            % disp([dd mm mean(newBounds)-mean(abs(dataTable_.bound(Ls)))])
            dataTable.bound(Ls) = sign(dataTable.bound(Ls)).*newBounds;
        end
    end
end

