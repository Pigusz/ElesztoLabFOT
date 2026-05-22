// --- Select input/output folders ---
input = getDirectory("Select the main folder containing the TXT files");
output = getDirectory("Select the output folder");

summaryPath = output + "summary_from_txt.csv";

// --- Preparation ---
if (File.exists(summaryPath)) File.delete(summaryPath);
File.append("Filename,TRITC,FITC,Overlap,TRITC_Threshold,FITC_Threshold\n", summaryPath);

// --- Start recursive processing ---
processFolder(input);

// ================================
// --- Recursive folder traversal ---
// ================================
function processFolder(dir) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        item = list[i];
        path = dir + item;
        if (File.isDirectory(path)) {
            processFolder(path);
        } else {
            if (endsWith(item, "_FITC_Intensity.txt")) {
                base = replace(item, "_FITC_Intensity.txt", "");
                fitcPath = dir + item;
                tritcPath = dir + base + "_TRITC_Intensity.txt";
                if (!File.exists(tritcPath)) {
                    print("⚠ Missing TRITC file: " + tritcPath);
                    continue;
                }

                // --- Automatic thresholds (Otsu) ---
                fitcThreshold = autoThreshold(fitcPath);
                tritcThreshold = autoThreshold(tritcPath);
                print("FITC threshold: " + fitcThreshold + " TRITC threshold: " + tritcThreshold);

                // --- Single channel counts ---
                fitcCount = countPixelsAboveThreshold(fitcPath, fitcThreshold);
                tritcCount = countPixelsAboveThreshold(tritcPath, tritcThreshold);

                // --- Simple overlap (both > threshold at same position) ---
                overlap = countOverlapPixels(fitcPath, tritcPath, fitcThreshold, tritcThreshold);

                // --- Save to CSV ---
                File.append(base + "," + tritcCount + "," + fitcCount + "," + overlap + "," + tritcThreshold + "," + fitcThreshold + "\n", summaryPath);
                print("✅ " + base);
            }
        }
    }
}

// ================================
// --- Otsu auto threshold ---
// ================================
function autoThreshold(path) {
    contents = File.openAsString(path);
    contents = replace(contents, "\r", "");
    lines = split(contents, "\n");

    // First pass: min, max, total pixels
    minVal = 1e300;
    maxVal = -1e300;
    total = 0;
    for (i = 0; i < lines.length; i++) {
        line = trim(lines[i]);
        if (line == "") continue;
        tokens = split(line, "\t");
        for (j = 0; j < tokens.length; j++) {
            v = parseFloat(trim(tokens[j]));
            if (!isNaN(v)) {
                if (v < minVal) minVal = v;
                if (v > maxVal) maxVal = v;
                total++;
            }
        }
    }

    if (total == 0 || maxVal == minVal) return maxVal;

    nBins = 256;
    binWidth = (maxVal - minVal) / nBins;
    histogram = newArray(nBins);
    for (i = 0; i < nBins; i++) histogram[i] = 0;

    // Fill histogram
    for (i = 0; i < lines.length; i++) {
        line = trim(lines[i]);
        if (line == "") continue;
        tokens = split(line, "\t");
        for (j = 0; j < tokens.length; j++) {
            v = parseFloat(trim(tokens[j]));
            if (!isNaN(v)) {
                bin = floor((v - minVal) / binWidth);
                if (bin >= nBins) bin = nBins - 1;
                histogram[bin]++;
            }
        }
    }

    // Otsu's method
    sumTotal = 0;
    for (i = 0; i < nBins; i++) sumTotal += i * histogram[i];
    sumB = 0; wB = 0; maxVariance = 0; thresholdBin = 0;
    for (t = 0; t < nBins; t++) {
        wB += histogram[t];
        if (wB == 0) continue;
        wF = total - wB;
        if (wF == 0) break;
        sumB += t * histogram[t];
        mB = sumB / wB;
        mF = (sumTotal - sumB) / wF;
        varBetween = wB * wF * (mB - mF) * (mB - mF);
        if (varBetween > maxVariance) {
            maxVariance = varBetween;
            thresholdBin = t;
        }
    }
    threshold = minVal + (thresholdBin + 1) * binWidth;
    return threshold;
}

// ================================
// --- Count pixels above threshold in a single file ---
// ================================
function countPixelsAboveThreshold(txtPath, threshold) {
    contents = File.openAsString(txtPath);
    contents = replace(contents, "\r", "");
    lines = split(contents, "\n");
    count = 0;
    for (j = 0; j < lines.length; j++) {
        line = trim(lines[j]);
        if (line == "") continue;
        values = split(line, "\t");
        for (k = 0; k < values.length; k++) {
            val = parseFloat(trim(values[k]));
            if (!isNaN(val) && val > threshold) count++;
        }
    }
    return count;
}

function countOverlapPixels(fitcPath, tritcPath, fitcThr, tritcThr) {
    // 1. Fájlok beolvasása
    fitcRaw = File.openAsString(fitcPath);
    tritcRaw = File.openAsString(tritcPath);

    // 2. Sorokra bontás
    fitcLines = split(replace(fitcRaw, "\r", ""), "\n");
    tritcLines = split(replace(tritcRaw, "\r", ""), "\n");

    // Diagnosztika
    print("Fájl beolvasva: " + fitcLines.length + " sor.");

    count = 0;
    maxF = 0;
    maxT = 0;

    // 3. Sorok bejárása
    for (i = 0; i < fitcLines.length; i++) {
        valsF = split(fitcLines[i], "\t");
        valsT = split(tritcLines[i], "\t");

        if (valsF.length < 2) valsF = split(fitcLines[i], " ");
        if (valsT.length < 2) valsT = split(tritcLines[i], " ");

        for (j = 0; j < valsF.length; j++) {

            vF = parseFloat(valsF[j]);
            vT = parseFloat(valsT[j]);


            if (!isNaN(vF) && vF > maxF) maxF = vF;
            if (!isNaN(vT) && vT > maxT) maxT = vT;

           
            if (!isNaN(vF) && !isNaN(vT)) {
                if (vF > fitcThr && vT > tritcThr) {
                    count++;
                }
            }
        }
    }

    // debug célokra
    print("--- Elemzés eredménye ---");
    print("Max érték (FITC): " + maxF);
    print("Max érték (TRITC): " + maxT);
    print("Talált átfedő pixelek: " + count);
    print("-------------------------");

    return count;
}

// ================================
// --- Utility min/max/length ---
// ================================
function minOf(a, b) { return (a < b) ? a : b; }
function maxOf(a, b) { return (a > b) ? a : b; }
function lengthOf(arr) { return arr.length; }