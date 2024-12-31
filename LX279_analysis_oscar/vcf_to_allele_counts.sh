#!/bin/bash

# Input and output file names
input_file="matrix_SNV.vcf"  # Initial VCF file
simplified_file="simplified.vcf"
filtered_file="filtered.vcf"
allele_reads_file="allele.reads.per.cell.tsv"
vaf_file="vaf.per.cell.tsv"
ref_output_file="snvs_scnano_ref.mtx"
alt_output_file="snvs_scnano_alt.mtx"

# Step 1: Simplify the VCF file

echo "--------------------------------------------------"
echo "Step 1. Simplifying VCF by keeping SNV position, base changed, calling QUAL and SNV counts per cell."

awk '
    BEGIN { FS = OFS = "\t" }
    /^#CHROM/ {
        sub(/^#/, "", $1);
        print;
        next;
    }
    /^#/ { next; }
    {
        for (i = 1; i <= NF; i++) {
            if (i != 7 && i != 8 && i != 9) {
                printf "%s%s", $i, (i == NF || (i == 6 && NF == 9) ? "\n" : "\t");
            }
        }
    }
' "$input_file" > "$simplified_file"
echo "Simplified VCF file created: $simplified_file"

echo "--------------------------------------------------"

# Step 2: Filter the VCF file

echo "Step 2. Removing multiallelic SNVs."
grep 'CHROM' "$simplified_file" > "$filtered_file"
awk '$5 ~ /^[ACTG]$/ {print}' "$simplified_file" >> "$filtered_file"
echo "Biallelic only VCF file created: $filtered_file"
echo "--------------------------------------------------"

# Step 3: Process the filtered VCF to generate allele reads per cell

echo "Step 3. Obtaining number of allele reads per cell."

awk '
    BEGIN { FS = OFS = "\t" }
    NR == 1 {
        print;
        next;
    }
    {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t", $1, $2, $3, $4, $5, $6;
        for (i = 7; i <= NF; i++) {
            split($i, fields, ":");
            printf "%s\t", fields[8];
        }
        print "";
    }
' "$filtered_file" > "$allele_reads_file"
echo "Allele reads per cell file created: $allele_reads_file"
echo "--------------------------------------------------"

# Step 4: Calculate VAF per cell

echo "Step 4. Calculating VAFs from allele reads per cell. NAs might be introduced."

awk '
    BEGIN { FS = OFS = "\t" }
    function process_genotype(genotype) {
        if (genotype == "0/0" || genotype == "0/0/0/0") return "NA";
        if (genotype ~ /0\/[0-9]+\/0\/0|0\/0\/[0-9]+\/0|0\/0\/0\/[0-9]+/) return "1";
        if (genotype ~ /^[0-9]+\/[0-9]+$/) {
            split(genotype, parts, "/");
            return parts[2] / (parts[1] + parts[2]);
        }
        if (genotype == "0/n") return "1";
        if (genotype == "n/0") return "0";
        if (genotype ~ /^[0-9]+\/[0-9]+\/[0-9]+\/[0-9]+$/) {
            split(genotype, parts, "/");
            sum = parts[2] + parts[3] + parts[4];
            return sum / (parts[1] + parts[2] + parts[3] + parts[4]);
        }
        return genotype;
    }
    NR == 1 {
        print;
        next;
    }
    {
        for (i = 6; i <= NF; i++) {
            $i = process_genotype($i);
        }
        print;
    }
' "$allele_reads_file" > "$vaf_file"
echo "VAF per cell file created: $vaf_file"
echo "--------------------------------------------------"

# Step 5: Divide the allele reads into reference and alternate files

echo "Step 5. Splitting allele reads per cell into ref/alt reads per cell."

awk '
BEGIN { FS = OFS = "\t" }
NR == 1 {
    print $0 > "'"$ref_output_file"'";
    print $0 > "'"$alt_output_file"'";
    next;
}
{
    ref_result = "";
    alt_result = "";
    for (i = 7; i <= NF; i++) {
        split($i, reads, "/");
        ref_result = (ref_result == "" ? reads[1] : ref_result OFS reads[1]);
        alt_result = (alt_result == "" ? reads[2] : alt_result OFS reads[2]);
    }
    ref_out = $1;
    alt_out = $1;
    for (i = 2; i <= 6; i++) {
        ref_out = ref_out OFS $i;
        alt_out = alt_out OFS $i;
    }
    print ref_out, ref_result > "'"$ref_output_file"'";
    print alt_out, alt_result > "'"$alt_output_file"'";
}
' "$allele_reads_file"
echo "Reference and alternate allele matrix files created: $ref_output_file and $alt_output_file"
echo "--------------------------------------------------"
