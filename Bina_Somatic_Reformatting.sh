#!/bin/bash

JAVA_1_7_DIR="/isilon/sequencing/Kurt/Programs/Java/jdk1.7.0_25/bin"
GATK_DIR="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.1-1"
REF_GENOME="/isilon/sequencing/GATK_resource_bundle/bwa_mem_0.7.5a_ref/human_g1k_v37_decoy.fasta"

BINA_DIR=$1 # UP TO THE HASH VALUE
OUTPUT_DIR=$2
SAMPLE_ID_PAIR=$3 # WHATEVER SAMPLE ID YOU WANT TO GIVE FOR THE PAIR
SAMPLE_ID_TUMOR=$4 # WHATEVER SAMPLE ID YOU WANT TO GIVE TO THE TUMOR

# $JAVA_1_7_DIR/java -jar \
# $GATK_DIR/GenomeAnalysisTK.jar \
# -T CombineVariants \
# -R $REF_GENOME \
# --variant $BINA_DIR/MergeSimple/InDel.Merge/variants.vcf \
# --variant $BINA_DIR/MergeSimple/SNP.Merge/variants.vcf \
# -o $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf"

( grep "^##" $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf" ; \
echo "##INFO=<ID=DP,Number=1,Type=String,Description=\"Dummy value so that annovar can be ran on this file\">" ; \
grep "^#CHROM" $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf" ; \
grep -v "^#" $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf" \
| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8";DP=100",$9,$10,$11}' ) \
>| $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.mod.vcf"

# ( grep "^##" $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf" ; \
# echo ##INFO=\<ID=DP,Number=1,Type=String,Description=\"Dummy value so that annovar can be ran on this file\"\> ; \
# grep "^#CHROM" $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf" ; \
# grep -v "^#" $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf" \
# | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8";DP=100",$9,$10,$11}' ) \
# >| $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.mod.vcf"

$JAVA_1_7_DIR/java -jar \
$GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--variant $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.mod.vcf" \
-sn TUMOR \
-o $OUTPUT_DIR/$SAMPLE_ID_TUMOR".bina.somatic.MergeSimple.combined.mod.vcf"
