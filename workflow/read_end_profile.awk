# FPR by distance to the nearest read end and by read-length quintile (supplement S1).
# USAGE: zcat calls.tsv.gz | awk -v OUT=<prefix> -v Q1=663 -v Q2=748 -v Q3=967 -v Q4=1296 -f read_end_profile.awk
function klass(k, b) { return (b == "A") ? "A" : ((substr(k, 4, 1) == "G") ? "CpG" : "CpH") }
function distbin(d) {
    if (d < 50)   return "1_0-50"
    if (d < 100)  return "2_50-100"
    if (d < 250)  return "3_100-250"
    if (d < 500)  return "4_250-500"
    if (d < 1000) return "5_500-1000"
    return "6_1000+"
}
function lenbin(L) {
    if (L <= Q1) return "1_<=" Q1
    if (L <= Q2) return "2_<=" Q2
    if (L <= Q3) return "3_<=" Q3
    if (L <= Q4) return "4_<=" Q4
    return "5_>" Q4
}
BEGIN {
    FS = OFS = "\t"
    nT = split("0.5,0.7,0.9,0.95,0.98,0.99", T, ",")
    if (OUT == "" || Q1 == "") { print "need OUT and Q1..Q4" > "/dev/stderr"; exit 2 }
}
NR == 1 { next }
{
    b = $18
    if (b != "C" && b != "A") next
    k = toupper($17)
    if (length(k) != 5 || k ~ /[^ACGT]/) { drop++; next }
    cls = klass(k, b)
    p = ($20 == "true") ? 0 : $13 + 0

    fp = $2 + 0; L = $12 + 0                      # forward_read_position, read_length
    d = fp; if (L - 1 - fp < d) d = L - 1 - fp
    db = distbin(d); lb = lenbin(L)

    dD[cls, db]++; dL[cls, lb]++
    seenD[db] = 1; seenL[lb] = 1
    for (i = 1; i <= nT; i++) if (p > T[i]) { nD[cls, db, i]++; nL[cls, lb, i]++ }
    tot++
}
END {
    fd = OUT ".readend.tsv"; fl = OUT ".readlen.tsv"
    hdr = "class" OFS "bin" OFS "n_obs"
    for (i = 1; i <= nT; i++) hdr = hdr OFS "FPR_" T[i]
    print hdr > fd; print hdr > fl
    # Bins are written in a fixed order.
    nd = split("1_0-50,2_50-100,3_100-250,4_250-500,5_500-1000,6_1000+", DB, ",")
    nl = split("1_<=" Q1 ",2_<=" Q2 ",3_<=" Q3 ",4_<=" Q4 ",5_>" Q4, LB, ",")
    split("CpG CpH A", ord, " ")
    for (o = 1; o <= 3; o++) {
        c = ord[o]
        for (j = 1; j <= nd; j++) { bn = DB[j]; if (dD[c, bn] == "") continue
            line = c OFS bn OFS dD[c, bn]
            for (i = 1; i <= nT; i++) line = line OFS sprintf("%.4e", (nD[c, bn, i] + 0) / dD[c, bn])
            print line > fd }
        for (j = 1; j <= nl; j++) { bn = LB[j]; if (dL[c, bn] == "") continue
            line = c OFS bn OFS dL[c, bn]
            for (i = 1; i <= nT; i++) line = line OFS sprintf("%.4e", (nL[c, bn, i] + 0) / dL[c, bn])
            print line > fl }
    }
    printf "rows=%d dropped=%d\n", tot, drop > "/dev/stderr"
}
