version 1.0

# =============================================================================
# MAC filter + rsID annotation
# =============================================================================
#
# One scatter shard per summary-statistics file. Each shard makes a single
# streaming pass that:
#   - merge-joins the sumstats with a variant annotation TSV on (chrom, pos),
#     matching ref/alt exactly, and prepends the annotation's rsid (NA if none)
#   - keeps only variants with MAF >= maf_threshold and
#     MAC = MAF * sample_size >= mac_threshold (rows with missing AF dropped)
# The output is bgzipped and tabix-indexed (.tbi).
#
# Both files must be sorted by chromosome (1..22, X/23, Y/24, MT/25), then
# position. The shard fails if either file goes out of order.
#
# Designed for FinnGen Refinery (Cromwell on GCP).
# =============================================================================

workflow mac_filter_rsidify {

  input {
    File    sumstats_manifest   # Headerless TSV: trait \t analysis_type \t gs://path \t sample_size
    File    annotation          # Variant annotation TSV (.gz ok) with chr, pos, ref, alt, rsid columns

    Int     mac_threshold = 30
    Float   maf_threshold = 0.0001

    String  docker

    # Summary-statistics column names
    String  chrom_col = "#chrom"
    String  pos_col   = "pos"
    String  ref_col   = "ref"
    String  alt_col   = "alt"
    String  af_col    = "af_alt"
  }

  scatter (row in read_tsv(sumstats_manifest)) {
    call rsidify_mac_filter {
      input:
        trait         = row[0],
        analysis_type = row[1],
        sumstats      = row[2],
        sample_size   = row[3],
        annotation    = annotation,
        mac_threshold = mac_threshold,
        maf_threshold = maf_threshold,
        chrom_col     = chrom_col,
        pos_col       = pos_col,
        ref_col       = ref_col,
        alt_col       = alt_col,
        af_col        = af_col,
        docker        = docker
    }
  }

  call gather_summary {
    input:
      stats  = rsidify_mac_filter.stats,
      docker = docker
  }

  output {
    Array[File] filtered_outputs = rsidify_mac_filter.filtered
    Array[File] filtered_tbi     = rsidify_mac_filter.filtered_tbi
    File        summary_tsv      = gather_summary.summary
  }
}


# =============================================================================
# Task: rsidify_mac_filter
# =============================================================================

task rsidify_mac_filter {

  input {
    String trait
    String analysis_type
    File   sumstats
    String sample_size
    File   annotation
    Int    mac_threshold
    Float  maf_threshold
    String chrom_col
    String pos_col
    String ref_col
    String alt_col
    String af_col
    String docker
  }

  String out_name = "~{trait}_~{analysis_type}_mac~{mac_threshold}_maf~{maf_threshold}_filtered.gz"

  command <<<
    set -euo pipefail

    gzip -dc -f "~{sumstats}" \
    | awk -F'\t' -v OFS='\t' \
        -v annot_cmd="gzip -dc -f '~{annotation}'" \
        -v trait="~{trait}" \
        -v analysis_type="~{analysis_type}" \
        -v sample_size="~{sample_size}" \
        -v mac_threshold="~{mac_threshold}" \
        -v maf_threshold="~{maf_threshold}" \
        -v chrom_col="~{chrom_col}" \
        -v pos_col="~{pos_col}" \
        -v ref_col="~{ref_col}" \
        -v alt_col="~{alt_col}" \
        -v af_col="~{af_col}" \
    '
    function die(msg) { print "ERROR: " msg > "/dev/stderr"; failed = 1; exit 1 }

    # Chromosome sort rank: 1..22, X=23, Y=24, MT=25
    function rank(c) {
      sub(/^chr/, "", c)
      if (c == "X") return 23
      if (c == "Y") return 24
      if (c == "M" || c == "MT") return 25
      return c + 0
    }

    function find(names, n, name,   i) {
      for (i = 1; i <= n; i++) if (names[i] == name) return i
      die("column \"" name "\" not found")
    }

    # Advance the annotation stream by one row
    function next_annot(   line, f) {
      if ((annot_cmd | getline line) <= 0) { a_eof = 1; return }
      split(line, f, "\t")
      a_r = rank(f[A_CHR]); a_p = f[A_POS] + 0
      if (a_r < a_prev_r || (a_r == a_prev_r && a_p < a_prev_p))
        die("annotation not sorted at " f[A_CHR] ":" f[A_POS])
      a_prev_r = a_r; a_prev_p = a_p
      a_key = f[A_REF] ":" f[A_ALT]; a_id = f[A_RSID]
    }

    BEGIN {
      if ((annot_cmd | getline line) <= 0) die("annotation file is empty")
      n = split(line, h, "\t")
      A_CHR = find(h, n, "chr"); A_POS = find(h, n, "pos")
      A_REF = find(h, n, "ref"); A_ALT = find(h, n, "alt"); A_RSID = find(h, n, "rsid")
      next_annot()
    }

    NR == 1 {
      n = split($0, h, "\t")
      S_CHR = find(h, n, chrom_col); S_POS = find(h, n, pos_col)
      S_REF = find(h, n, ref_col);   S_ALT = find(h, n, alt_col); S_AF = find(h, n, af_col)
      # Output columns are shifted by one (rsid first): tell tabix where chrom/pos are
      print S_CHR + 1, S_POS + 1 > "tabix_cols.txt"
      sub(/^#/, "", $S_CHR)   # "#chrom" -> "chrom"; "#rsid" is now the header marker
      print "#rsid", $0
      next
    }

    {
      total++
      r = rank($S_CHR); p = $S_POS + 0
      if (r < s_prev_r || (r == s_prev_r && p < s_prev_p))
        die("sumstats not sorted at " $S_CHR ":" $S_POS)
      s_prev_r = r; s_prev_p = p

      # New position: skip annotation rows before it, buffer alleles at it
      if (r != buf_r || p != buf_p) {
        delete buf; buf_r = r; buf_p = p
        while (!a_eof && (a_r < r || (a_r == r && a_p < p))) next_annot()
        while (!a_eof && a_r == r && a_p == p) { buf[a_key] = a_id; next_annot() }
      }

      # MAC/MAF filter
      af = $S_AF
      if (af == "" || af == "NA") next
      af += 0
      maf = (af > 0.5) ? 1 - af : af
      if (maf < maf_threshold || maf * sample_size < mac_threshold) next

      key = $S_REF ":" $S_ALT
      id = (key in buf) ? buf[key] : "NA"
      if (id != "NA") with_rsid++
      kept++
      print id, $0
    }

    END {
      if (failed) exit 1
      close(annot_cmd)
      printf "%s\t%s\t%s\t%d\t%d\t%d\n", trait, analysis_type, sample_size, total, kept, with_rsid > "stats.tsv"
      printf "%s: %d variants, %d pass MAC/MAF, %d of those with rsid\n", trait, total, kept, with_rsid > "/dev/stderr"
    }' \
    | bgzip -@ 2 -c > "~{out_name}"

    # Header line starts with "#rsid", so tabix treats it as the header
    read -r seq_col pos_col < tabix_cols.txt
    tabix -s "${seq_col}" -b "${pos_col}" -e "${pos_col}" "~{out_name}"
  >>>

  output {
    File   filtered     = out_name
    File   filtered_tbi = out_name + ".tbi"
    String stats        = read_string("stats.tsv")
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "2 GB"
    disks: "local-disk " + (ceil(size(sumstats, "GB") * 3 + size(annotation, "GB")) + 5) + " HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
    noAddress: true
  }
}


# =============================================================================
# Task: gather_summary — one row per input file
# =============================================================================

task gather_summary {

  input {
    Array[String] stats
    String        docker
  }

  command <<<
    set -euo pipefail
    printf 'trait\tanalysis_type\tsample_size\ttotal_variants\tfiltered_variants\twith_rsid\n' > summary.tsv
    cat ~{write_lines(stats)} >> summary.tsv
  >>>

  output {
    File summary = "summary.tsv"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
    noAddress: true
  }
}
