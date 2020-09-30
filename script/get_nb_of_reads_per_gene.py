import pysam
import matplotlib.pyplot as plt

bam = "BYK_GDB_ONT_1_PAD67469_Aflipflop.against_transcriptome.secondary_and_primary_only.sorted.bam"
genes_dict = {}
for ali in pysam.AlignmentFile(bam, 'rb'):
    if not ali.is_secondary :
        gene_name = ali.reference_name.split("|")[5]
        transcript_name = ali.reference_name.split("|")[4]
        if gene_name not in genes_dict.keys():
            genes_dict[gene_name] = {transcript_name : 1}
        else :
            if transcript_name in genes_dict[gene_name].keys():
                genes_dict[gene_name][transcript_name] += 1
            else:
                genes_dict[gene_name][transcript_name] = 1

print(genes_dict["Tmod2"])

groups = [0,1,10,20,50,100]

nb_reads_per_genes = {}
for gene, transcripts in genes_dict.items():
    nb_reads = sum([count for count in transcripts.values()])
    nb_reads_per_genes[gene] = nb_reads