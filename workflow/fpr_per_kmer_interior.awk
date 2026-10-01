# Per-5-mer false-positive table from modkit 0.6.4 `extract full`, read interiors (>= 250 bp from a read end).
# USAGE: zcat calls.tsv.gz | awk -v OUT=<prefix> -f fpr_per_kmer_interior.awk  ->  <prefix>.{ref,query}.tsv
# Inferred rows count as observations with p = 0; ref k-mers are case-normalised and strand-corrected. See docs/methods_notes.md.
function revcomp(s,   i, c, o) {
    o = ""
    for (i = length(s); i >= 1; i--) {
        c = substr(s, i, 1)
        o = o ((c == "A") ? "T" : (c == "C") ? "G" : (c == "G") ? "C" : (c == "T") ? "A" : "N")
    }
    return o
}
function klass(k, b,   nxt) {
    if (b == "A") return "A"
    nxt = substr(k, 4, 1)
    return (nxt == "G") ? "CpG" : "CpH"
}
BEGIN {
    FS = OFS = "\t"
    nT = split("0.5,0.6,0.7,0.8,0.9,0.95,0.98,0.99", T, ",")
    if (OUT == "") { print "set -v OUT=<prefix>" > "/dev/stderr"; exit 2 }
}
NR == 1 { next }
{
    # Read-end exclusion: d = min(fw_pos, read_length - 1 - fw_pos); rows with d < 250 leave numerator and denominator.
    fwp = $2 + 0; rlen = $12 + 0
    dend = (fwp < (rlen - 1 - fwp)) ? fwp : (rlen - 1 - fwp)
    if (dend < 250) { zoneskip++; next }

    b = $18                                       # canonical base
    if (b != "C" && b != "A") { odd++; next }
    p = ($20 == "true") ? 0 : $13 + 0             # inferred -> probability 0

    rk = toupper($16)                             # soft-masked reference -> normalise case
    if ($7 == "-") rk = revcomp(rk)
    qk = toupper($17)

    if (length(rk) != 5 || rk ~ /[^ACGT]/) { dropR++; skip_r = 1 } else skip_r = 0
    if (length(qk) != 5 || qk ~ /[^ACGT]/) { dropQ++; skip_q = 1 } else skip_q = 0
    if (skip_r && skip_q) next

    if (substr(rk, 3, 1) == b) refctr_ok++         # audit counters
    if (substr(qk, 3, 1) == b) qryctr_ok++
    if (rk != qk) mismatch++
    n++

    if (!skip_r) { dR[rk]++; if (!(rk in clsR)) clsR[rk] = klass(rk, b) }
    if (!skip_q) { dQ[qk]++; if (!(qk in clsQ)) clsQ[qk] = klass(qk, b) }
    for (i = 1; i <= nT; i++) {
        if (p > T[i]) { if (!skip_r) nR[rk, i]++; if (!skip_q) nQ[qk, i]++ }
    }
}
END {
    fr = OUT ".ref.tsv"; fq = OUT ".query.tsv"
    hdr = "kmer" OFS "class" OFS "n_obs"
    for (i = 1; i <= nT; i++) hdr = hdr OFS "fp_" T[i] OFS "FPR_" T[i] OFS "Q_" T[i]
    print hdr > fr
    print hdr > fq

    for (k in dR) emit(k, dR[k], clsR[k], "R", fr)
    for (k in dQ) emit(k, dQ[k], clsQ[k], "Q", fq)

    printf "rows_used=%d  odd_base_rows=%d  dropped_outer250=%d\n", n, odd, zoneskip                        > "/dev/stderr"
    printf "dropped (non-ACGT / short kmer): ref %d (%.3f%%), query %d (%.3f%%)\n", \
           dropR, 100*dropR/(n+dropR), dropQ, 100*dropQ/(n+dropQ)            > "/dev/stderr"
    printf "distinct kmers: ref %d, query %d  (max possible 512)\n", length(dR), length(dQ) > "/dev/stderr"
    printf "ref_kmer centre == canonical after strand fix: %.2f%%\n", 100*refctr_ok/n > "/dev/stderr"
    printf "query_kmer centre == canonical:                %.2f%%\n", 100*qryctr_ok/n > "/dev/stderr"
    printf "ref_kmer != query_kmer after strand fix:       %.2f%%\n", 100*mismatch/n  > "/dev/stderr"
}
# n < 200 -> NA. A context with no false positive gets Q from a single-event pseudocount, marked '>' as a bound.
function emit(k, nobs, cls, which, f,   i, fp, fpr, q, line) {
    line = k OFS cls OFS nobs
    for (i = 1; i <= nT; i++) {
        fp = (which == "R") ? nR[k, i] + 0 : nQ[k, i] + 0
        if (nobs < 200) { line = line OFS fp OFS "NA" OFS "NA"; continue }
        fpr = fp / nobs
        if (fp == 0) { q = -10 * log(1 / nobs) / log(10); line = line OFS fp OFS "0" OFS sprintf(">%.1f", q) }
        else         { q = -10 * log(fpr)     / log(10); line = line OFS fp OFS sprintf("%.3e", fpr) OFS sprintf("%.1f", q) }
    }
    print line > f
}
