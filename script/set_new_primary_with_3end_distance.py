import pysam
import sys
from random import choice

## A bam file sorted by read name is needed : uncomment pysam.sort line if bam file isn't sorted yet.
    # E.g : pysam.sort("-n -o", sortedbam, bamfile)

def pick_best_AS_candidates(alignments):
    ## Return list of alignments with the best alignment score
    best_AS = max([alignment.get_tag("AS") for alignment in alignments])
    best_candidates = [alignment for alignment in alignments if alignment.get_tag("AS") == best_AS]
    return best_candidates

def get_dist_from_3end(alignment, sortedbam):
    ## Return the distance (in number of base) between the 3' end of the transcript and the 3' end of the read.
    transcript_ID = alignment.reference_name
    transcript_length = sortedbam.get_reference_length(transcript_ID)
    assert transcript_length >= alignment.reference_end, "Negative length !"
    dist_from_3end = transcript_length - alignment.reference_end
    # print("transcript " + alignment.reference_name + " : " + alignment.query_name + " : " +str(transcript_dict[alignment.reference_name]) + " | " + str(alignment.reference_end))
    return dist_from_3end

def pick_best_3end_dist(alignments, sortedbam):
    ## Return list of alignment with the shortest distance between the 3' end of the transcript and the 3' end of the read.
    best_candidates = []
    dist_from_3end_list = [get_dist_from_3end(alignment, sortedbam) for alignment in alignments]
    min_dist_from_3end = min(dist_from_3end_list)
    for i, dist_from_3end in enumerate(dist_from_3end_list):
        if dist_from_3end == min_dist_from_3end:
            best_candidates.append(alignments[i])
    return best_candidates

def pick_best_candidate(alignments, sortedbam):
    ## Return one alignment with max alignment score, shortest distance to 3' end of the transcript, then a random one in case of equalities.
    best_candidates = pick_best_AS_candidates(alignments)
    if len(best_candidates) > 1 :
        best_candidates = pick_best_3end_dist(best_candidates, sortedbam)
    best_candidate = choice(best_candidates)
    return best_candidate

def set_flag_to_not_supplementary(alignment):
    alignment.flag = int(bin(alignment.flag & 2047),2)

def set_flag_to_primary(alignment):
    alignment.flag = int(bin(alignment.flag & 3839),2)

def set_flag_to_secondary(alignment):
    alignment.flag = int(bin(alignment.flag | 256),2)

def set_new_flag(alignments, best_alignment):
    ## Set the best alignment flag (in regard of AS and shortest distance to 3' end of the transcript) to 'primary', secondary for the others.
    for alignment in alignments :
        set_flag_to_not_supplementary(alignment)
        if alignment == best_alignment :
            set_flag_to_primary(alignment)
        else :
            set_flag_to_secondary(alignment)

def output_new_primary_alignments(alignments, output, sortedbam):
    ## For every alignments of one read :
    ##     1) set the primary alignment flag to the alignment with best AS and shortest distance between its 3' end and the 3' end of the transcript.
    ##     2) set the secondary flag to the others
    ##     3) write them in a new bam file
    if alignments :
        # if len(alignments) > 1 :
        best_candidate = pick_best_candidate(alignments, sortedbam)
            # if best_candidate.is_secondary :
        set_new_flag(alignments, best_candidate)
            # print("New best alignment for read " + best_candidate.query_name + " : " + str(best_candidate.flag))
        for alignment in alignments :
            if not alignment.is_secondary:
                output.write(alignment)

def has_QC_greater_than_75(alignment):
    lengthRead = alignment.infer_query_length()
    lengthAlignment = alignment.query_alignment_length
    percentAlignment = float(lengthAlignment + 20)/ float(lengthRead)
    return percentAlignment>=0.80

def main() :
    inputFile = sys.argv[1]
    if len(sys.argv) > 2 :
        outputFile = sys.argv[2]
    else :
        outputFileName = inputFile[:-3] + "filter_80QC_new_primary.bam"
    # sorted_inputFile = inputFile[:-3] + "sorted_by_read_name.bam"
    # pysam.sort("-n", inputFile, "-o", sorted_inputFile)
    # inputFile = sorted_inputFile
    inputFile = pysam.AlignmentFile(inputFile, 'rb')
    outputFile = pysam.AlignmentFile(outputFileName, "wb", template = inputFile)
    read_name = None
    read_alignments = []
    for alignment in inputFile :
        if not alignment.is_unmapped and not alignment.is_supplementary and has_QC_greater_than_75(alignment):
            if alignment.query_name != read_name :
                output_new_primary_alignments(read_alignments, outputFile, inputFile)
                read_alignments = []
                read_name = alignment.query_name
            read_alignments.append(alignment)
        else :
            pass
    output_new_primary_alignments(read_alignments, outputFile, inputFile) # Catch the last read

if __name__ == "__main__":
    main()
