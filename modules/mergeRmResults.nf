/*
========================================================================================
    MERGE REPEATMASKER RESULTS
========================================================================================
    Combine adjusted RepeatMasker results from all batches into final outputs
    that match the original REPEAT_MASK module interface.

    Steps:
      1. Concatenate all adjusted .out files, renumber linkage IDs globally,
         and recompute (left) column using authoritative scaffold lengths
      2. Reassemble masked FASTA — merge chunks back into full scaffolds in order
      3. Generate GFF from merged .out
      4. Generate standard RM-format summary .tbl from the .out file
      5. Generate detailed by-family .tbl using buildSummary.pl
      6. Merge .align files, compute divergence, generate repeat landscape

    Inputs:
        - assembly_id:     Sample identifier
        - out_files:       Collection of adjusted .out files
        - fasta_files:     Collection of adjusted masked FASTA files
        - align_files:     Collection of adjusted .align files
        - seq_lengths:     Scaffold lengths TSV (from PREP_RM_BATCHES)

    Outputs:
        - masked_fasta:    Soft-masked genome
        - repeat_gff:      Repeat annotation GFF
        - repeat_tbl:      Standard RM-format summary table
        - repeat_out:      Merged RepeatMasker .out annotation file
        - repeat_tbl_full: Detailed by-family table (buildSummary.pl)
        - repeat_align:    Merged .align file
        - repeat_divsum:   Divergence summary from calcDivergenceFromAlign.pl
        - repeat_landscape: Repeat landscape HTML plot

    Tools: RepeatMasker utilities (buildSummary.pl, rmOutToGFF3.pl,
           calcDivergenceFromAlign.pl, createRepeatLandscape.pl)
           from tetools container
========================================================================================
*/

process MERGE_RM_RESULTS {
    tag "${assembly_id}"
    label 'merge_rm'
    publishDir "${params.outdir}/repeats", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(out_files), path(fasta_files), path(align_files)
    tuple val(assembly_id), path(seq_lengths)

    output:
    tuple val(assembly_id), path("${assembly_id}.masked.fasta"),         emit: masked_fasta
    tuple val(assembly_id), path("${assembly_id}.repeat.gff3"),                  emit: repeat_gff
    tuple val(assembly_id), path("${assembly_id}.tbl"),                  emit: repeat_tbl
    tuple val(assembly_id), path("${assembly_id}.out"),                  emit: repeat_out
    tuple val(assembly_id), path("${assembly_id}.buildSummary.tbl"),     emit: repeat_tbl_full
    tuple val(assembly_id), path("${assembly_id}.align"),                emit: repeat_align
    tuple val(assembly_id), path("${assembly_id}.divsum"),               emit: repeat_divsum
    tuple val(assembly_id), path("${assembly_id}.repeat_landscape.html"), emit: repeat_landscape

    script:
    """
    set -euo pipefail

    # ── 1. Merge .out files with globally renumbered IDs ──────────────────────
    head -3 \$(ls -1 ${out_files} | sort -V | head -1) > merged_raw.out

    ls -1 ${out_files} | sort -V | while read -r f; do
        tail -n +4 "\$f" >> merged_raw.out
    done

    # Renumber linkage IDs (column 15) and recompute (left) column (column 8)
    # using authoritative scaffold lengths to eliminate cross-batch inconsistencies
    awk -v lenfile="${seq_lengths}" '
    BEGIN {
        while ((getline line < lenfile) > 0) {
            split(line, a, "\\t")
            slen[a[1]] = int(a[2])
        }
        close(lenfile)
    }
    NR <= 3 { print; next }
    {
        sid = \$5
        if (sid in slen) {
            \$8 = "(" slen[sid] - \$7 ")"
        }
        \$15 = ++counter
        print
    }
    ' merged_raw.out > ${assembly_id}.out

    # ── 2. Reassemble masked FASTA ───────────────────────────────────────────
    > chunk_index.tsv
    for fa in ${fasta_files}; do
        grep "^>" "\$fa" | while read -r hdr; do
            name=\${hdr#>}
            scaffold=\$(echo "\$name" | sed 's/__[0-9]*__[0-9]*\$//')
            start=\$(echo "\$name" | grep -oP '__\\K[0-9]+(?=__[0-9]+\$)')
            echo -e "\${scaffold}\\t\${start}\\t\${fa}\\t\${name}"
        done
    done >> chunk_index.tsv

    sort -k1,1 -k2,2n chunk_index.tsv > chunk_index_sorted.tsv

    > ${assembly_id}.masked.fasta
    current_scaffold=""
    while IFS=\$'\\t' read -r scaffold start fa seqname; do
        if [ "\$scaffold" != "\$current_scaffold" ]; then
            if [ -n "\$current_scaffold" ]; then
                echo "" >> ${assembly_id}.masked.fasta
            fi
            echo ">\${scaffold}" >> ${assembly_id}.masked.fasta
            current_scaffold="\$scaffold"
        fi
        awk -v target=">\${seqname}" '
            BEGIN { found=0 }
            /^>/ { found = (\$0 == target) ? 1 : 0; next }
            found { printf "%s", \$0 }
        ' "\$fa" >> ${assembly_id}.masked.fasta
    done < chunk_index_sorted.tsv
    echo "" >> ${assembly_id}.masked.fasta

    # Reformat to standard FASTA line width
    awk '
    /^>/ { if (NR > 1) print ""; print; next }
    { printf "%s", \$0 }
    END { print "" }
    ' ${assembly_id}.masked.fasta \\
    | fold -w 60 > tmp_wrapped.fa
    mv tmp_wrapped.fa ${assembly_id}.masked.fasta

    # ── 3. Generate GFF3 from .out with proper SO feature types ──────────────
    python3 - ${assembly_id}.out ${assembly_id}.repeat.gff3 <<'PYEOF'
import sys, re

out_file, gff_file = sys.argv[1], sys.argv[2]

# RepeatMasker class/family → Sequence Ontology feature type
SO_MAP = {
    "SINE":             "SINE_element",
    "LINE":             "LINE_element",
    "LTR":              "LTR_retrotransposon",
    "DNA":              "DNA_transposon",
    "RC":               "helitron",
    "Retroposon":       "retroposon",
    "Simple_repeat":    "simple_repeat",
    "Low_complexity":   "low_complexity_region",
    "Satellite":        "satellite_DNA",
    "rRNA":             "rRNA",
    "tRNA":             "tRNA",
    "snRNA":            "snRNA",
    "srpRNA":           "SRP_RNA",
    "scRNA":            "scRNA",
    "Unknown":          "repeat_region",
}

def classify(class_family):
    # Map RM class/family string to SO term
    main = class_family.split("/")[0]
    return SO_MAP.get(main, "repeat_region")

with open(out_file) as fin, open(gff_file, "w") as fout:
    fout.write("##gff-version 3\\n")
    for i, line in enumerate(fin):
        if i < 3:
            continue
        fields = line.split()
        if len(fields) < 15:
            continue

        sw_score    = fields[0]
        perc_div    = fields[1]
        perc_del    = fields[2]
        perc_ins    = fields[3]
        scaffold    = fields[4]
        qstart      = int(fields[5])
        qend        = int(fields[6])
        strand_raw  = fields[8]
        repeat_name = fields[9]
        class_fam   = fields[10]
        hit_id      = fields[14]

        strand = "-" if strand_raw == "C" else "+"
        so_type = classify(class_fam)

        # Target coordinates: always forward-orientation (start < end)
        # + strand .out cols: repeat_begin(11) repeat_end(12) (repeat_left)(13)
        # C strand .out cols: (repeat_left)(11) repeat_end(12) repeat_begin(13)
        if strand == "+":
            target_start = fields[11]
            target_end   = fields[12]
        else:
            target_start = fields[13]
            target_end   = fields[12]

        attrs = (
            f"ID=repeat_{hit_id};"
            f"Name={repeat_name};"
            f"class={class_fam};"
            f"pctdiv={perc_div};"
            f"pctdel={perc_del};"
            f"pctins={perc_ins};"
            f"Target={repeat_name} {target_start} {target_end} {strand}"
        )

        fout.write(f"{scaffold}\\tRepeatMasker\\t{so_type}\\t{qstart}\\t{qend}\\t{sw_score}\\t{strand}\\t.\\t{attrs}\\n")

PYEOF

    # ── 4. Generate standard RM-format .tbl from .out ────────────────────────
    perl -e '
    use strict; use warnings;

    my (\$total_len, \$n_excl, \$n_seqs, \$gc_count, \$bases_masked) = (0,0,0,0,0);
    open my \$fa, "<", "\$ARGV[0]" or die;
    while (<\$fa>) {
        chomp;
        if (/^>/) { \$n_seqs++; next }
        \$total_len += length(\$_);
        my \$s = uc(\$_);
        \$gc_count += (\$s =~ tr/GC//);
        my \$acgt = \$_;
        \$acgt =~ s/[^ACGTacgt]//g;
        \$n_excl += length(\$acgt);
        my \$lc = \$_;
        \$lc =~ s/[^acgt]//g;
        \$bases_masked += length(\$lc);
    }
    close \$fa;

    my \$gc_pct = \$n_excl > 0 ? (\$gc_count/\$n_excl)*100 : 0;
    my \$mask_pct = \$total_len > 0 ? (\$bases_masked/\$total_len)*100 : 0;

    my (%main_count, %main_bp, %sub_count, %sub_bp);
    my (\$interspersed_n, \$interspersed_bp) = (0, 0);

    open my \$out, "<", "\$ARGV[1]" or die;
    my \$ln = 0;
    while (<\$out>) {
        \$ln++;
        next if \$ln <= 3;
        chomp;
        my @f = split;
        next unless @f >= 11;
        my \$class_fam = \$f[10];
        my \$bp = \$f[6] - \$f[5] + 1;

        my (\$main, \$sub) = split /\\//, \$class_fam, 2;
        \$sub //= "";

        \$main_count{\$main} += 1;
        \$main_bp{\$main} += \$bp;
        if (\$sub ne "") {
            \$sub_count{"\$main/\$sub"} += 1;
            \$sub_bp{"\$main/\$sub"} += \$bp;
        }

        my %non_interspersed = map { \$_ => 1 } qw(Simple_repeat Low_complexity Satellite rRNA tRNA snRNA srpRNA scRNA);
        unless (\$non_interspersed{\$main}) {
            \$interspersed_n++;
            \$interspersed_bp += \$bp;
        }
    }
    close \$out;

    my \$pct = sub { \$total_len > 0 ? (\$_[0]/\$total_len)*100 : 0 };

    printf "==================================================\\n";
    printf "file name: merged\\n";
    printf "sequences:         %5d\\n", \$n_seqs;
    printf "total length: %10d bp  (%d bp excl N/X-runs)\\n", \$total_len, \$n_excl;
    printf "GC level:       %6.2f %%\\n", \$gc_pct;
    printf "bases masked: %10d bp ( %5.2f %%)\\n", \$bases_masked, \$mask_pct;
    printf "==================================================\\n";
    printf "               number of      length   percentage\\n";
    printf "               elements*    occupied  of sequence\\n";
    printf "--------------------------------------------------\\n";

    my \$retro_n  = (\$main_count{LTR}//0) + (\$main_count{LINE}//0) + (\$main_count{SINE}//0);
    my \$retro_bp = (\$main_bp{LTR}//0)    + (\$main_bp{LINE}//0)    + (\$main_bp{SINE}//0);
    printf "Retroelements     %5d   %10d bp  %5.2f %%\\n", \$retro_n, \$retro_bp, \$pct->(\$retro_bp);
    printf "   SINEs:         %5d   %10d bp  %5.2f %%\\n", \$main_count{SINE}//0, \$main_bp{SINE}//0, \$pct->(\$main_bp{SINE}//0);
    printf "   LINEs:         %5d   %10d bp  %5.2f %%\\n", \$main_count{LINE}//0, \$main_bp{LINE}//0, \$pct->(\$main_bp{LINE}//0);
    printf "   LTR elements:  %5d   %10d bp  %5.2f %%\\n", \$main_count{LTR}//0, \$main_bp{LTR}//0, \$pct->(\$main_bp{LTR}//0);

    for my \$k (sort keys %sub_count) {
        my (\$m, \$s) = split /\\//, \$k, 2;
        next unless \$m eq "LTR";
        printf "     %-12s %5d   %10d bp  %5.2f %%\\n", \$s, \$sub_count{\$k}, \$sub_bp{\$k}, \$pct->(\$sub_bp{\$k});
    }

    printf "\\n";
    printf "DNA transposons   %5d   %10d bp  %5.2f %%\\n", \$main_count{DNA}//0, \$main_bp{DNA}//0, \$pct->(\$main_bp{DNA}//0);
    for my \$k (sort keys %sub_count) {
        my (\$m, \$s) = split /\\//, \$k, 2;
        next unless \$m eq "DNA";
        printf "   %-14s %5d   %10d bp  %5.2f %%\\n", \$s, \$sub_count{\$k}, \$sub_bp{\$k}, \$pct->(\$sub_bp{\$k});
    }

    printf "\\n";
    printf "Unclassified:     %5d   %10d bp  %5.2f %%\\n", \$main_count{Unknown}//0, \$main_bp{Unknown}//0, \$pct->(\$main_bp{Unknown}//0);
    printf "\\nTotal interspersed repeats:  %10d bp  %5.2f %%\\n\\n", \$interspersed_bp, \$pct->(\$interspersed_bp);

    my \$sm_n  = (\$main_count{rRNA}//0)+(\$main_count{tRNA}//0)+(\$main_count{snRNA}//0)+(\$main_count{srpRNA}//0)+(\$main_count{scRNA}//0);
    my \$sm_bp = (\$main_bp{rRNA}//0)+(\$main_bp{tRNA}//0)+(\$main_bp{snRNA}//0)+(\$main_bp{srpRNA}//0)+(\$main_bp{scRNA}//0);
    printf "Small RNA:        %5d   %10d bp  %5.2f %%\\n", \$sm_n, \$sm_bp, \$pct->(\$sm_bp);
    printf "Satellites:       %5d   %10d bp  %5.2f %%\\n", \$main_count{Satellite}//0, \$main_bp{Satellite}//0, \$pct->(\$main_bp{Satellite}//0);
    printf "Simple repeats:   %5d   %10d bp  %5.2f %%\\n", \$main_count{Simple_repeat}//0, \$main_bp{Simple_repeat}//0, \$pct->(\$main_bp{Simple_repeat}//0);
    printf "Low complexity:   %5d   %10d bp  %5.2f %%\\n", \$main_count{Low_complexity}//0, \$main_bp{Low_complexity}//0, \$pct->(\$main_bp{Low_complexity}//0);
    printf "==================================================\\n";
    ' ${assembly_id}.masked.fasta ${assembly_id}.out > ${assembly_id}.tbl

    # ── 5. Generate detailed by-family table (buildSummary.pl) ────────────────
    #   Use the authoritative seq_lengths (from seqkit fx2tab) directly as the
    #   -genome parameter. This avoids ACGT-only counting from the masked FASTA
    #   which excludes Ns and causes length mismatches with the .out file.
    perl /opt/RepeatMasker/util/buildSummary.pl \
        -genome ${seq_lengths} \
        ${assembly_id}.out > ${assembly_id}.buildSummary.tbl

    # ── 6. Merge .align files and generate repeat landscape ───────────────────
    cat \$(ls -1 ${align_files} | sort -V) > ${assembly_id}.align

    # Calculate Kimura 2-parameter divergence
    perl /opt/RepeatMasker/util/calcDivergenceFromAlign.pl \
        -s ${assembly_id}.divsum \
        -a ${assembly_id}.align_with_div \
        ${assembly_id}.align

    # Generate repeat landscape HTML
    genome_size=\$(awk '!/^>/{total += length(\$0)} END{print total}' ${assembly_id}.masked.fasta)

    perl /opt/RepeatMasker/util/createRepeatLandscape.pl \
        -div ${assembly_id}.divsum \
        -g \${genome_size} \
        > ${assembly_id}.repeat_landscape.html
    """

    stub:
    """
    touch ${assembly_id}.masked.fasta
    touch ${assembly_id}.repeat.gff3
    touch ${assembly_id}.tbl
    touch ${assembly_id}.out
    touch ${assembly_id}.buildSummary.tbl
    touch ${assembly_id}.align
    touch ${assembly_id}.divsum
    touch ${assembly_id}.repeat_landscape.html
    """
}