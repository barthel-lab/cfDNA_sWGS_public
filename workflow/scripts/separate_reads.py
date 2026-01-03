#!/usr/bin/env python

import sys
import re
import pysam

def is_mouse_chrom(chrom):
    if chrom is None:
        return False
    return bool(re.search(r'(_mm|\.)', chrom))

def is_human_chrom(chrom):
    if chrom is None:
        return False
    return not is_mouse_chrom(chrom)

def check_ambiguous(xa_tag, primary_chrom):
    # Returns True if XA tag indicates cross-species alignment
    if is_human_chrom(primary_chrom) and re.search(r'chr.*_mm', xa_tag):
        return True
    if is_mouse_chrom(primary_chrom) and re.search(r'chr[^_]*[^.]($|[^_])', xa_tag):
        return True
    return False

def separate_reads(input_bam_path, output_human_bam_path, output_mouse_bam_path, output_amb_bam_path):
    bam_in = pysam.AlignmentFile(input_bam_path, "rb")
    
    # Open outputs with the same header as input BAM
    bam_out_human = pysam.AlignmentFile(output_human_bam_path, "wb", template=bam_in)
    bam_out_mouse = pysam.AlignmentFile(output_mouse_bam_path, "wb", template=bam_in)
    bam_out_amb = pysam.AlignmentFile(output_amb_bam_path, "wb", template=bam_in)
    
    #!/usr/bin/env python

import sys
import re
import pysam

def is_mouse_chrom(chrom):
    if chrom is None:
        return False
    return bool(re.search(r'(_mm|\.)', chrom))

def is_human_chrom(chrom):
    if chrom is None:
        return False
    return not is_mouse_chrom(chrom)

def check_ambiguous(xa_tag, sa_tag, primary_chrom):
    # Check cross-species hits in XA or SA tags
    # XA: alternative hits from primary
    if xa_tag:
        if is_human_chrom(primary_chrom) and re.search(r'chr.*_mm', xa_tag):
            return True
        if is_mouse_chrom(primary_chrom) and re.search(r'chr[^_]*[^.]($|[^_])', xa_tag):
            return True

    # SA: supplementary alignments
    if sa_tag:
        if is_human_chrom(primary_chrom) and re.search(r'chr.*_mm', sa_tag):
            return True
        if is_mouse_chrom(primary_chrom) and re.search(r'chr[^_]*[^.]($|[^_])', sa_tag):
            return True

    return False

def separate_reads(input_bam_path, output_human_bam_path, output_mouse_bam_path, output_amb_bam_path):
    bam_in = pysam.AlignmentFile(input_bam_path, "rb")
    
    # Open outputs with the same header as input BAM
    bam_out_human = pysam.AlignmentFile(output_human_bam_path, "wb", template=bam_in)
    bam_out_mouse = pysam.AlignmentFile(output_mouse_bam_path, "wb", template=bam_in)
    bam_out_amb = pysam.AlignmentFile(output_amb_bam_path, "wb", template=bam_in)
    
    for read in bam_in.fetch(until_eof=True):
        primary_chrom = read.reference_name
        mate_chrom = None
        if not read.is_unmapped and not read.mate_is_unmapped:
            mate_chrom = bam_in.get_reference_name(read.next_reference_id)
        
        # If primary read unmapped, classify ambiguous
        if primary_chrom is None:
            bam_out_amb.write(read)
            continue
        
        # Check XA and SA tags for cross-species alignment
        xa_tag = None
        sa_tag = None
        try:
            xa_tag = read.get_tag("XA")
        except KeyError:
            pass
        try:
            sa_tag = read.get_tag("SA")
        except KeyError:
            pass

        if check_ambiguous(xa_tag, sa_tag, primary_chrom):
            bam_out_amb.write(read)
            continue

        
        # Exclude secondary/supplementary alignments
        if read.is_secondary or read.is_supplementary:
            # Optional: decide if you want to ignore or assign elsewhere
            continue
        
        # CONFIDENT mouse: primary & mate on mouse chr or mate unmapped
        if is_mouse_chrom(primary_chrom) and (mate_chrom is None or is_mouse_chrom(mate_chrom)):
            bam_out_mouse.write(read)
            continue
        
        # CONFIDENT human: primary & mate on human chr or mate unmapped
        if is_human_chrom(primary_chrom) and (mate_chrom is None or is_human_chrom(mate_chrom)):
            bam_out_human.write(read)
            continue
        
        # Otherwise ambiguous
        bam_out_amb.write(read)
    
    bam_in.close()
    bam_out_human.close()
    bam_out_mouse.close()
    bam_out_amb.close()

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: separate_reads.py input.bam output_human.bam output_mouse.bam output_amb.bam")
        sys.exit(1)
    
    separate_reads(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
