#!/usr/bin/env python3

import argparse
import pandas as pd
from collections import defaultdict
from typing import Dict, List, Set, Tuple, NamedTuple

class SequenceInfo(NamedTuple):
    """Container for sequence classification information"""
    is_t2t: bool
    chrom_name: str
    level: str
    sequence_length: int
    num_gaps: int
    total_gap_length: int

def read_mashmap(mashmap_file: str) -> pd.DataFrame:
    """
    Read and process mashmap output file.
    """
    cols = ['query', 'query_length', 'query_start', 'query_end', 'strand',
            'target', 'target_length', 'target_start', 'target_end', 'matches',
            'block_len', 'qual', 'identity', 'kmer_score']
    return pd.read_csv(mashmap_file, sep='\t', names=cols)

def read_bed(bed_file: str) -> Dict[str, List[Tuple[int, int]]]:
    """Read BED file and return dictionary of intervals."""
    intervals = defaultdict(list)
    with open(bed_file, 'r') as f:
        for line in f:
            chrom, start, end = line.strip().split()[:3]
            intervals[chrom].append((int(start), int(end)))
    return intervals

def is_alignment_valid(alignment: pd.Series, 
                      ignored_regions: Dict[str, List[Tuple[int, int]]]) -> bool:
    """
    Check if alignment is valid (has parts outside ignored regions).
    Returns True if either:
    1. Alignment is not on a chromosome with ignored regions
    2. Alignment extends beyond ignored regions
    """
    if alignment['target'] not in ignored_regions:
        return True
    
    align_start = alignment['target_start']
    align_end = alignment['target_end']
    
    # Check if alignment extends beyond any ignored region
    for start, end in ignored_regions[alignment['target']]:
        if align_start < start or align_end > end:
            return True
    return False

def has_telomeres(seq_id: str, seq_length: int, 
                  telomere_regions: Dict[str, List[Tuple[int, int]]]) -> bool:
    """
    Check if sequence has telomeres at both ends.
    A sequence is considered to have telomeres if there are telomeric repeats
    within 10kb of both ends.
    """
    if seq_id not in telomere_regions:
        return False
    
    tel_regions = sorted(telomere_regions[seq_id])
    if len(tel_regions) < 2:
        return False
    
    return tel_regions[0][0] < 10000 and tel_regions[-1][1] > seq_length - 10000

def process_sequence(seq_id: str,
                    seq_length: int,
                    seq_aligns: pd.DataFrame,
                    ignored_regions: Dict[str, List[Tuple[int, int]]],
                    telomere_regions: Dict[str, List[Tuple[int, int]]],
                    gap_regions: Dict[str, List[Tuple[int, int]]]) -> SequenceInfo:
    """
    Process a single sequence and determine its classification.
    
    The classification logic considers:
    1. Valid alignments (outside ignored regions)
    2. Coverage and uniqueness of mapping
    3. Presence of telomeres for T2T determination
    4. Gap content for scaffold/contig designation
    """
    # Filter to valid alignments
    valid_aligns = seq_aligns[seq_aligns.apply(
        lambda x: is_alignment_valid(x, ignored_regions), axis=1)]
    
    # Calculate gap statistics
    num_gaps = len(gap_regions.get(seq_id, []))
    total_gap_length = sum(end - start for start, end in gap_regions.get(seq_id, []))
    level = "scaffold" if num_gaps > 0 else "contig"
    
    # Extract Genbank ID from sequence ID (format: sample_id#hap_int#genbank_id)
    genbank_id = seq_id.split('#')[2]
    
    # If no valid alignments exist, classify as chrUn
    if len(valid_aligns) == 0:
        return SequenceInfo(
            is_t2t=False,
            chrom_name=f"chrUn_{genbank_id}",
            level=level,
            sequence_length=seq_length,
            num_gaps=num_gaps,
            total_gap_length=total_gap_length
        )
    
    # Calculate coverage for each target chromosome using block_len
    coverages = {}
    for target in valid_aligns['target'].unique():
        target_aligns = valid_aligns[valid_aligns['target'] == target]
        total_covered = target_aligns['block_len'].sum()
        coverages[target] = total_covered / seq_length
    
    # Find primary chromosome (one with highest coverage)
    primary_chr = max(coverages.items(), key=lambda x: x[1])[0]
    primary_coverage = coverages[primary_chr]
    
    # Check if it's T2T
    is_t2t = (primary_coverage > 0.95 and 
              has_telomeres(seq_id, seq_length, telomere_regions) and
              len(valid_aligns['target'].unique()) == 1)
    
    if is_t2t:
        chrom_name = primary_chr
    elif len(valid_aligns['target'].unique()) == 1:
        chrom_name = f"{primary_chr}_{genbank_id}_random"
    else:
        chrom_name = f"chrUn_{genbank_id}"
    
    return SequenceInfo(
        is_t2t=is_t2t,
        chrom_name=chrom_name,
        level=level,
        sequence_length=seq_length,
        num_gaps=num_gaps,
        total_gap_length=total_gap_length
    )

def process_alignments(mashmap_df: pd.DataFrame, 
                      ignored_regions: Dict[str, List[Tuple[int, int]]],
                      telomere_regions: Dict[str, List[Tuple[int, int]]],
                      gap_regions: Dict[str, List[Tuple[int, int]]]) -> Tuple[Dict[str, str], List[Dict]]:
    """
    Process all alignments to determine chromosome assignments and T2T status.
    Returns two items:
    1. Dictionary mapping sequence IDs to their chromosome assignments
    2. List of T2T sequences with their detailed information
    """
    # Get all sequence lengths from the first occurrence in mashmap data
    sequence_lengths = dict(mashmap_df.groupby('query')['query_length'].first())
    
    chrom_assignments = {}
    t2t_sequences = []
    
    # Process each sequence
    for seq_id, seq_length in sequence_lengths.items():
        # Get all alignments for this sequence
        seq_aligns = mashmap_df[mashmap_df['query'] == seq_id]
        
        # Process the sequence
        seq_info = process_sequence(
            seq_id, seq_length, seq_aligns,
            ignored_regions, telomere_regions, gap_regions
        )
        
        # Store chromosome assignment for all sequences
        chrom_assignments[seq_id] = seq_info.chrom_name
        
        # Only store T2T sequences in the t2t_sequences list
        if seq_info.is_t2t:
            t2t_sequences.append({
                'sequence_id': seq_id,
                'chr_name': seq_info.chrom_name,
                'level': seq_info.level,
                'sequence_length': seq_info.sequence_length,
                'num_gaps': seq_info.num_gaps,
                'total_gap_length': seq_info.total_gap_length
            })
    
    return chrom_assignments, t2t_sequences

def main():
    parser = argparse.ArgumentParser(description='Process assembly alignments and create annotation files')
    parser.add_argument('--mashmap', required=True, help='Mashmap alignment file')
    parser.add_argument('--telomeres', required=True, help='Telomere regions BED file')
    parser.add_argument('--gaps', required=True, help='Gap regions BED file')
    parser.add_argument('--ignore', required=True, help='Ignored regions BED file')
    parser.add_argument('--out-prefix', required=True, help='Output files prefix')
    
    args = parser.parse_args()
    
    # Read input files
    mashmap_df = read_mashmap(args.mashmap)
    telomere_regions = read_bed(args.telomeres)
    gap_regions = read_bed(args.gaps)
    ignored_regions = read_bed(args.ignore)
    
    # Process alignments
    chrom_assignments, t2t_sequences = process_alignments(
        mashmap_df, ignored_regions, telomere_regions, gap_regions)
    
    # Write chromAlias file with header and Genbank IDs
    with open(f"{args.out_prefix}.chromAlias.txt", 'w') as f:
        # Write header
        f.write("# assembly\tucsc\tgenbank\n")
        
        # Write data lines with extracted Genbank ID
        for seq_id, chrom_name in sorted(chrom_assignments.items()):
            # Extract Genbank ID from sequence ID (format: sample_id#hap_int#genbank_id)
            genbank_id = seq_id.split('#')[2]
            f.write(f"{seq_id}\t{chrom_name}\t{genbank_id}\n")
    
    # Write T2T chromosomes file (only for T2T sequences)
    if t2t_sequences:
        pd.DataFrame(t2t_sequences).to_csv(
            f"{args.out_prefix}.t2t_chromosomes.tsv", 
            sep='\t', 
            index=False
        )

if __name__ == "__main__":
    main()