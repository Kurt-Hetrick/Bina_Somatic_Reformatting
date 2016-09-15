#!/bin/bash



JAVA_1_7_DIR="/isilon/sequencing/Kurt/Programs/Java/jdk1.7.0_25/bin"

# JAVA_1_7_DIR="/isilon/sequencing/Kurt/Programs/Java/jdk1.8.0_73/bin"

GATK_DIR="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.1-1"

# GATK_DIR="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.6"
REF_GENOME="/isilon/sequencing/GATK_resource_bundle/bwa_mem_0.7.5a_ref/human_g1k_v37_decoy.fasta"

BINA_DIR=$1 # UP TO THE HASH VALUE
OUTPUT_DIR=$2
SAMPLE_ID_PAIR=$3 # WHATEVER SAMPLE ID YOU WANT TO GIVE FOR THE PAIR (VCF FILE NAME)
SAMPLE_ID_TUMOR=$4 # WHATEVER SAMPLE ID YOU WANT TO GIVE TO THE TUMOR (VCF FILE NAME)
SAMPLE_ID_NORMAL=$5 # WHATEVER SAMPLE ID YOU WANT TO GIVE THE NORMAL (VCF FILE NAME)
TUMOR_GERMLINE_VCF=$6 # THE FULL VCF FILE PATH FOR THE REGULAR VCF FILE
NORMAL_GERMLINE_VCF=$7 # THE FULL VCF FILE PATH FOR THE REGULAR VCF FILE

# Merge the InDel adn SNP outputs from Bina consensus caller

$JAVA_1_7_DIR/java -jar \
$GATK_DIR/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $REF_GENOME \
--variant $BINA_DIR/MergeSimple/InDel.Merge/variants.vcf \
--variant $BINA_DIR/MergeSimple/SNP.Merge/variants.vcf \
-o $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf"

# Add a dummy INFO header for Depth
# Add a dummy value for DP in the INFO tag
# Needed so that annovar can be ran

( grep "^##" $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf" ; \
echo "##INFO=<ID=DP,Number=1,Type=String,Description=\"Dummy value so that annovar can be ran on this file\">" ; \
grep "^#CHROM" $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf" ; \
grep -v "^#" $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.vcf" \
| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8";DP=100",$9,$10,$11}' ) \
>| $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.mod.vcf"

# Grab the tumor out of the paired vcf

$JAVA_1_7_DIR/java -jar \
$GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--variant $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.mod.vcf" \
-sn TUMOR \
-o $OUTPUT_DIR/$SAMPLE_ID_TUMOR".bina.somatic.MergeSimple.combined.mod.vcf"

# Grab the normal out of the paired vcf

$JAVA_1_7_DIR/java -jar \
$GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
--variant $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.mod.vcf" \
-sn NORMAL \
-o $OUTPUT_DIR/$SAMPLE_ID_NORMAL".bina.somatic.MergeSimple.combined.mod.vcf"

# TAKE THE PAIRED VCF AND TURN IT INTO A TAB DELIMTED TABLE

$JAVA_1_7_DIR/java -jar \
$GATK_DIR/GenomeAnalysisTK.jar \
-T VariantsToTable \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--variant $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.mod.vcf" \
--fields CHROM \
--fields POS \
--fields ID \
--fields REF \
--fields ALT \
--fields QUAL \
--fields EVENTLENGTH \
--fields FILTER \
--fields HOM-REF \
--fields HET \
--fields HOM-VAR \
--fields NO-CALL \
--fields SOMATIC \
--fields SOURCES \
--genotypeFields AD \
--genotypeFields GT \
--genotypeFields DP \
--genotypeFields DP4 \
--genotypeFields RD \
--genotypeFields SC \
--genotypeFields SHIFT3 \
--allowMissingData \
--showFiltered \
-o $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.mod.txt"

# CLASSIFY INDELS AND SNVS
# CONVERT THE NUCLEOTIDE LEVEL GENOTYPES INTO MORE GENERALIZED HOM_REF, HET AND HOM_VAR CLASSES

( head -n 1 $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.mod.txt" \
| awk 'BEGIN {OFS="\t"} {print $0,"VARIANT_TYPE","NORMAL_GT","TUMOR_GT"}' ; \
awk 'NR>1' $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.mod.txt" \
| awk '$7==0 {print $0"\t""SNV"} $7!="0" {print $0"\t""INDEL"}' \
| awk '{split($16,NORMAL_GT,"/"); if (NORMAL_GT[1]==NORMAL_GT[2]&&NORMAL_GT[1]==$4) print $0"\t""HOM_REF"; \
else if (NORMAL_GT[1]!=NORMAL_GT[2]) print $0"\t""HET" ; \
else print $0"\t""HOM_VAR"}' \
| awk '{split($23,TUMOR_GT,"/"); if (TUMOR_GT[1]==TUMOR_GT[2]&&TUMOR_GT[1]==$4) print $0"\t""HOM_REF"; \
else if (TUMOR_GT[1]!=TUMOR_GT[2]) print $0"\t""HET" ; \
else print $0"\t""HOM_VAR"}') \
>| $OUTPUT_DIR/$SAMPLE_ID_PAIR".bina.somatic.MergeSimple.combined.table.txt"

# make a temp directory to create vcf files that have the sequence dictionary removed from the header.

mkdir -p $OUTPUT_DIR/TEMP/

( grep "^#" $OUTPUT_DIR/$SAMPLE_ID_NORMAL".bina.somatic.MergeSimple.combined.mod.vcf" \
| grep -v "^##contig=" ; \
grep -v "^#" $OUTPUT_DIR/$SAMPLE_ID_NORMAL".bina.somatic.MergeSimple.combined.mod.vcf" \
| egrep -v "^GL|^NC|^hs37d5") \
>| $OUTPUT_DIR/TEMP/$SAMPLE_ID_NORMAL".bina.somatic.MergeSimple.combined.mod.vcf"

# Duke's vcf files had a few tags in the INFO section that was not defined in the header. so had to define them.
# Also had to remove calls on the non-primary assembly contigs

(grep "^##" $NORMAL_GERMLINE_VCF \
| grep -v "^##contig=" ; \
echo "##INFO=<ID=AnnoDBID,Number=1,Type=Integer,Description=\"Something that Duke put in here\">" ; \
echo "##FILTER=<ID=QDFilter,Description=\"Undefined filter that Duke uses\">" ; \
echo "##FILTER=<ID=FSFilter,Description=\"Undefined filter that Duke uses\">" ; \
grep "^#CHROM" $NORMAL_GERMLINE_VCF \
| cut -f 1-9 \
| awk 'BEGIN {OFS="\t"} {print $0,"NORMAL"}' ; \
grep -v "^#" $NORMAL_GERMLINE_VCF \
| egrep -v "^GL|^NC|^hs37d5") \
>| $OUTPUT_DIR/TEMP/NORMAL_GERMLINE_VCF.vcf

( grep "^#" $OUTPUT_DIR/$SAMPLE_ID_TUMOR".bina.somatic.MergeSimple.combined.mod.vcf" \
| grep -v "^##contig=" ; \
grep -v "^#" $OUTPUT_DIR/$SAMPLE_ID_TUMOR".bina.somatic.MergeSimple.combined.mod.vcf" \
| egrep -v "^GL|^NC|^hs37d5") \
>| $OUTPUT_DIR/TEMP/$SAMPLE_ID_TUMOR".bina.somatic.MergeSimple.combined.mod.vcf"

# Duke's vcf files had a few tags in the INFO section that was not defined in the header. so had to define them.
# Also had to remove calls on the non-primary assembly contigs

(grep "^##" $TUMOR_GERMLINE_VCF \
| grep -v "^##contig=" ; \
echo "##INFO=<ID=AnnoDBID,Number=1,Type=Integer,Description=\"Something that Duke put in here\">" ; \
echo "##FILTER=<ID=QDFilter,Description=\"Undefined filter that Duke uses\">" ; \
echo "##FILTER=<ID=FSFilter,Description=\"Undefined filter that Duke uses\">" ; \
grep "^#CHROM" $TUMOR_GERMLINE_VCF \
| cut -f 1-9 \
| awk 'BEGIN {OFS="\t"} {print $0,"TUMOR"}' ; \
grep -v "^#" $TUMOR_GERMLINE_VCF \
| egrep -v "^GL|^NC|^hs37d5") \
>| $OUTPUT_DIR/TEMP/TUMOR_GERMLINE_VCF.vcf

# Merge the germline calls with the filtered calls for the normal sample
# If there is a simple collision, keep the germline call
# If there is complex collision, keep both

$JAVA_1_7_DIR/java -jar \
$GATK_DIR/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $REF_GENOME \
--variant:bina $OUTPUT_DIR/TEMP/$SAMPLE_ID_NORMAL".bina.somatic.MergeSimple.combined.mod.vcf" \
--variant:germline $OUTPUT_DIR/TEMP/NORMAL_GERMLINE_VCF.vcf \
--rod_priority_list germline,bina \
--genotypemergeoption PRIORITIZE \
--printComplexMerges \
-o $OUTPUT_DIR/$SAMPLE_ID_NORMAL".bina_plus_germline.vcf"

# Merge the germline calls with the filtered calls for the tumor sample
# If there is a simple collision, keep the germline call
# If there is complex collision, keep both

$JAVA_1_7_DIR/java -jar \
$GATK_DIR/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $REF_GENOME \
--variant:bina $OUTPUT_DIR/TEMP/$SAMPLE_ID_TUMOR".bina.somatic.MergeSimple.combined.mod.vcf" \
--variant:germline $OUTPUT_DIR/TEMP/TUMOR_GERMLINE_VCF.vcf \
--rod_priority_list germline,bina \
--genotypemergeoption PRIORITIZE \
--printComplexMerges \
-o $OUTPUT_DIR/$SAMPLE_ID_TUMOR".bina_plus_germline.vcf"
