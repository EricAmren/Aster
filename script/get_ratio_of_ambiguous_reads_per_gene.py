import HTSeq
import pysam
import csv

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# BAM FILES NEED TO BE SORTED BY READ NAME !!

# bam = "BYK_GDB_ONT_1_PAD67469_Aflipflop.against_transcriptome.secondary_and_primary_only.sorted.bam"
bam = "BYK_GDB_ONT_1_PAD67469_Aflipflop.against_transcriptome.secondary_and_primary_only.sorted_by_read_name.bam"
# bam = "BYK_GCC_ONT_1_FAK82381_A.basecallHA.primary_and_secondary_only.sorted_by_read_names.bam"


def get_gene_name(ali):
    return ali.reference_name.split("|")[5]

def get_transcript_name(ali):
    return ali.reference_name.split("|")[4]

def bundle_multiple_alignments(alignment_file) :
    alignment_iter = iter(alignment_file)
    ali = next(alignment_iter)
    multiple_alignments = [ali]
    for ali in alignment_iter:
        if ali.query_name != multiple_alignments[0].query_name:
            yield multiple_alignments
            multiple_alignments = [ali]
        else :
            multiple_alignments.append(ali)
    yield multiple_alignments


def check_read_for_ambiguous_alignments(alignments, gene_dict):
    best_AS = max([alignment.get_tag("AS") for alignment in alignments])
    best_alignments = [alignment for alignment in alignments if alignment.get_tag("AS") == best_AS]
    increment_number_of_alignments(gene_dict, best_alignments)
    return gene_dict

            
def increment_number_of_alignments(gene_dict, alignments):
    mapped_transcripts = sorted([get_transcript_name(alignment) for alignment in alignments])
    transcripts_key = ""
    for transcript_name in mapped_transcripts:
        transcripts_key += transcript_name + "|"
    transcripts_key = transcripts_key[:-1]
    # transcripts_key = str(mapped_transcripts).replace("'","")
    mapped_genes = set(get_gene_name(alignment) for alignment in alignments)
    for gene in mapped_genes:
        if gene not in gene_dict.keys():
            gene_dict[gene] = {transcripts_key : 1}
        else :
            if transcripts_key not in gene_dict[gene].keys():
                gene_dict[gene][transcripts_key] = 1
            else :
                gene_dict[gene][transcripts_key] += 1
    return gene_dict



def output_gene_dict(gene_dict, output_file, read_threshold):
    with open(output_file, "w") as output:
        # gene_dict = {k : v for k, v in sorted(gene_dict.items(), key=lambda item: item[1])}
        output.write("Gene_name\tratio_of_ambiguous_reads\tnumber_of_reads\ttranscripts_names\tnumber_of_reads (ratio of uniquely mapped reads)\tmixed_transcripts\n")
        for gene, transcript_dict in gene_dict.items():
            nb_of_reads = sum([v for v in transcript_dict.values()])
            nb_of_ambiguous_reads = 0
            if nb_of_reads >= read_threshold :
                transcripts_str = ""
                transcripts_mix_str = ""
                nb_of_ambiguous_reads_dict = {}

                for transcript_key, count in transcript_dict.items():
                    transcript_name_list = transcript_key.split("|")
                    if len(transcript_name_list) > 1:
                        transcripts_mix_str += transcript_key + " : " + str(count) + "; "
                        nb_of_ambiguous_reads += count
                        for transcript in transcript_name_list :
                            if transcript not in nb_of_ambiguous_reads_dict.keys():
                                nb_of_ambiguous_reads_dict[transcript] = count
                            else :
                                nb_of_ambiguous_reads_dict[transcript] += count
                is_mapped_uniquely = False
                for transcript_key, count in transcript_dict.items(): # Count here is nb of non_ambiguous 
                    transcript_name_list = transcript_key.split("|")
                    if len(transcript_name_list) == 1:
                        is_mapped_uniquely = True
                        transcript = transcript_name_list[0]
                        if transcript in nb_of_ambiguous_reads_dict.keys():
                            nb_of_ambiguous_reads = nb_of_ambiguous_reads_dict[transcript]
                        else :
                            nb_of_ambiguous_reads = 0
                        ratio_of_non_ambiguous_reads_per_transcript = str(round(1 - nb_of_ambiguous_reads / (count + nb_of_ambiguous_reads) ,2))
                        transcripts_str += transcript_name_list[0] + " : " + str(count) + " (" + ratio_of_non_ambiguous_reads_per_transcript + "); "

                ratio = str(round(float(nb_of_ambiguous_reads)/float(nb_of_reads), 2)) + "\t"
                if is_mapped_uniquely :
                    transcripts_str = transcripts_str[:-2] + "\t"
                else :
                    transcripts_str = "(No read uniquely support this gene)\t"
                transcripts_mix_str = transcripts_mix_str[:-2] + "\t"
                output.write("%s\t%s%s\t%s%s\n"%(gene, ratio, nb_of_reads, transcripts_str, transcripts_mix_str))



                    

                # if counts[2]:
                #     sorted_transcripts_dict = {k : v for k, v in sorted(counts[2].items(), key=lambda item: item[1])}
                #     for transcripts_group, trans_count in sorted_transcripts_dict.items():
                #         transcript_str += str(transcripts_group) + " : " + str(trans_count) + "\t"

                # transcripts = str(counts[2]) if counts[2] else ""
                # ratio = round(float(counts[1])/float(counts[0]), 2)
                # output.write("%s\t%s\t%s\t%s\n"%(gene, str(ratio), str(counts[0]), transcript_str))
                # plot_data.append(ratio)

def main():
    alignment_file = pysam.AlignmentFile(bam)
    # new_bam = bam[:-3] + "only_ambiguous_alignment.bam"
    # new_bam_file = pysam.AlignmentFile(new_bam, 'wb', template = alignment_file)
    gene_dict = {}
    n = 0
    for alignment_bundle in bundle_multiple_alignments(alignment_file):
        n += 1
        check_read_for_ambiguous_alignments(alignment_bundle, gene_dict)
        if n % 100000 == 0 : 
            print(str(n) + " reads processed...")
            break
    output_file = 'ratio_of_ambiguous_reads_per_gene_100reads_GDB_TEST.csv'
    print("Writing results...")
    output_gene_dict(gene_dict, output_file, 100)
    print("Finished")
    
    # plot_data = []

    # sns.set(color_codes=True)
    # sns.distplot(plot_data, kde = False)
    # plt.show()

main()

