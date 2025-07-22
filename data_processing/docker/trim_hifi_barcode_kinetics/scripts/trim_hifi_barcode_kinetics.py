"""
trim_hifi_barcode_kinetics.py

Trim PacBio kinetics arrays (fi/fp/ri/rp) to match adapter‑clipped HiFi reads.

This script may be necessary to run for files which were demultiplexed using lima <2.5
These files retained kinetics values for the barcodes but the sequence values did not
resulting in a mismatch in their array lengths (and importantly a shift in the correct)
position for the kinetics value in relation to the base prediction.

Input
-----
uBAM or BAM from Sequel II/IIe that has:
    • adapters already trimmed
       • per‑read barcode‑clip lengths in tag ``bx:B:i,left,right``
    • raw kinetics in tags ``fi fp ri rp`` (B:C arrays)

Output
------
A BAM with the kinetics tags trimmed so that
``len(fi) == len(fp) == len(ri) == len(rp) == query_length``.
Optionally removes existing ``MM/ML`` methylation tags 


Example
-------
$ python trim_hifi_barcode_kinetics.py \
    in_name.hifi_reads.bam \
    out_name.hifi_reads.bam \
    --drop_mm \
    --threads 8
"""

from __future__ import annotations

import argparse
import sys
import logging
import csv
from datetime import datetime
import array
from collections.abc import Sequence
import pysam
import os


class TrimError(RuntimeError):
    """Raised when a read cannot be trimmed consistently."""


def parse_bx(rec: pysam.AlignedSegment) -> tuple[int, int]:
    """
    Uses Pysam's get_tag() to pull the bx tag from a read. Note that the barcode tag
    is expected to be of the form "bx:B:i,left,right" where left and right are 
    probably 18 (or something similar).

    Parameters:
        rec (pysam.AlignedSegment): read from PacBio uBam (which was demultiplexed)
    
    Returns:
        (left_trim, right_trim): integer trim values for left and right of read
    """
    try:
        bx_vals, _vtype = rec.get_tag("bx", with_value_type=True)
    except KeyError:
        raise TrimError("bx tag missing")

    if len(bx_vals) == 2:
        left, right = bx_vals

    else:
        raise TrimError(f"bx tag has unexpected length {len(bx_vals)}")

    if left < 0 or right < 0:
        raise TrimError("bx trim counts must be non-negative")

    return left, right


def trim_array(values: array.array | list[int] | list[float],
               left: int, right: int) -> array.array | list[int] | list[float]:
    """Return a slice removing *left* from the front and *right* from the end."""
    if right == 0:
        trimmed = values[left:]
    else:
        trimmed = values[left: -right]
    return trimmed


def remove_barcode_kinetics(rec: pysam.AlignedSegment,
             drop_mm: bool=True,
             metrics: dict[str, int] = None,
             read_logger: csv.writer = None) -> pysam.AlignedSegment:
    """
    Removes kinetics values from barcodes that were retained by some lima versions

    Parameters:
        rec (pysam.AlignedSegment): read from PacBio uBam (which was demultiplexed)
        drop_mm (bool): drop MM/ML tags (if present)
        metrics (dict[str, int], optional ): Dictionary of processing counters. N.B. mutated in-place.
            Expected keys include:
                - 'reads_skipped': Number of reads that didn't require trimming.
                - 'reads_fixed': Number of reads trimmed successfully.
                - 'mm_ml_removed': Number of reads where MM/ML tags were removed.
        read_logger (csv.writer; optional): A CSV writer used to log per-read trimming diagnostics.
            If provided, writes one row per read with the read name, sequence length,
            and the before/after lengths of fi, fp, ri, and rp kinetics tags.

    Returns:
        (pysam.AlignedSegment): read with barcode kinetics stripped out of kinetics array.


    """
    ## check to see what barcode length is from Bx tag
    left_trim, right_trim = parse_bx(rec)
    sequence_length = rec.query_length

    forward_tags = ("fi", "fp")
    reverse_tags = ("ri", "rp")
    all_tags = forward_tags + reverse_tags
    
    # Initialize per-read metrics
    log_data = {
        "read_name": rec.query_name,
        "sequence_length": sequence_length,
    }
    for tag in all_tags:
        log_data[f"{tag}_bf"] = -1
        log_data[f"{tag}_after"] = -1

    def log_read():
        """Emit one row to read_logger if active."""
        if read_logger:
            read_logger.writerow([
                log_data["read_name"],
                log_data["sequence_length"],
                log_data["fi_bf"], log_data["fi_after"],
                log_data["ri_bf"], log_data["ri_after"],
                log_data["fp_bf"], log_data["fp_after"],
                log_data["rp_bf"], log_data["rp_after"],
            ])

    # Skip if no trimming is needed
    if left_trim == 0 and right_trim == 0:
        for tag in all_tags:
            try:
                values, _ = rec.get_tag(tag, with_value_type=True)
                tag_len = len(values)
                log_data[f"{tag}_bf"] = tag_len
                log_data[f"{tag}_after"] = tag_len
            except KeyError:
                continue
        log_read()
        metrics["reads_skipped"] += 1
        return rec
    

    # Process each tag
    for tag in all_tags:
        try:
            values, vtype = rec.get_tag(tag, with_value_type=True)
        except KeyError:
            raise TrimError(f"required kinetics tag {tag} missing")

        if not isinstance(values, (array.array, list)) or isinstance(values, (str, bytes)):
            raise TrimError(f"kinetics tag {tag} is not a numeric array-like sequence")

        before_len = len(values)
        log_data[f"{tag}_bf"] = before_len

        ## Note that reverse-orientation tags will be trimmed in reverse.
        # So if the barcode trim  was 17, 18 then we will trim ri and rp as 18, 17
        ltrim, rtrim = (
            (left_trim, right_trim) if tag in forward_tags
            else (right_trim, left_trim)
        )

        ## see if positions to trim are larger than values to trim (it happens, see below)
        if ltrim + rtrim >= before_len:
            if before_len == 0:
                ##  According to PacBio's docs, sometimes one orientation may get filtered out 
                ## and empty lists are present for that orientation
                ## See https://pacbiofileformats.readthedocs.io/en/11.0/BAM.html

                ## We also see instances where neither orientation has kinetics values.
                ## In both cases just keep the empty kinetics arrays as this is acceptable.                
                new_vals = values  # valid empty array
            else:
                raise TrimError(f"{tag} trim ({ltrim},{rtrim}) exceeds array length {before_len}")
        else:
            ## Trim the kinetics array
            new_vals = trim_array(values, ltrim, rtrim)

        after_len = len(new_vals)
        log_data[f"{tag}_after"] = after_len
        
        ## set the trimmed values for the read
        rec.set_tag(tag, new_vals)



    ## Check to see if any tags are both: nonzero and have array lengths that differ in size
    ## from the length of the read (which might imply double trimming, for instance)
    length_mismatch = [
        f"{tag}={log_data[f'{tag}_after']}"
        for tag in all_tags
        if log_data[f"{tag}_after"] not in (0, sequence_length)
    ]

    ## exit!
    if length_mismatch:
        raise TrimError(
            f"kinetics length mismatch after trimming "
            f"(read len {sequence_length}; {', '.join(length_mismatch)})"
        )
      
    # Remove MM/ML if requested
    if drop_mm:
        for tag in ("MM", "ML"):
            if rec.has_tag(tag):
                ## I don't know how to delete a tag, but if you set with 
                ## value of None the tag is deleted. 
                rec.set_tag(tag, None)

                ## note we will remove twice as many MM/ML tags as reads
                ## so the count is set to 0.5 to prevent double counting
                metrics["mm_ml_removed"] += 0.5

    metrics["reads_fixed"] += 1
    log_read()

    return rec


def update_bam_header(in_header: pysam.AlignmentHeader, cmdline: str) -> pysam.AlignmentHeader:
    """
    Adds documentation for this script as a @PG line to pre-existing header.  Note this does not 
    update the header in place, instead a AlignmentHeader is rebuilt from a dict and returned.

    Parameters:
        in_header (pysam.AlignmentHeader): header from input file (that is getting fixed)
        cmdline (str): command line to append in @PG line

    Returns:
        (pysam.AlignmentHeader): updated header with @PG line documenting this script
    """

    hdr_dict = in_header.to_dict()

    # Add our program line
    pg_record = {
        "ID": "fix_pacbio_kinetics",
        "PN": "fix_pacbio_kinetics",
        "VN": "1.0",
        "CL": cmdline,
        "PP": hdr_dict["PG"][-1]["ID"] if hdr_dict.get("PG") else None,
        "DS": f"Trimmed fi/fp/ri/rp by bx on {datetime.now().isoformat(timespec='seconds')}",
    }
    hdr_dict.setdefault("PG", []).append(pg_record)

    return pysam.AlignmentHeader.from_dict(hdr_dict)



def process_file(in_path: str, 
                 out_path: str,
                 threads: int,
                 drop_mm: bool,
                 read_log_path: str = None) -> None:

    """
    Processes a PacBio BAM file to trim barcode-induced kinetic artifacts
    (fi/fp/ri/rp tags) and optionally remove MM/ML tags. Writes output to out_path.

    Parameters:
        in_path (str): Path to the input HiFi uBAM file (with kinetics and barcode info)
        out_path (str): Path where the output BAM file will be written.
        threads (int): Number of threads to use for BAM reading/writing.
        drop_mm (bool): If True, remove MM/ML tags from each read.
        read_log_path (str, optional): path/name for CSV for per-read kinetics tag trimming info
    """
    
    metrics = {
        "reads_total": 0,
        "reads_fixed": 0,
        "reads_skipped": 0,
        "reads_failed": 0,
        "mm_ml_removed": 0,
    }

    # Capture exact command line for header provenance
    cmdline = " ".join(sys.argv)

    # Set up per-read logging if requested
    read_log_file = None
    read_logger = None
    if read_log_path:
        read_log_file = open(read_log_path, 'w', newline='')
        read_logger = csv.writer(read_log_file)
        # Write log file header
        read_logger.writerow([
            'read_name', 'sequence_length',
            'fi_bf', 'fi_after',
            'ri_bf', 'ri_after', 
            'fp_bf', 'fp_after',
            'rp_bf', 'rp_after'
        ])
        logging.info(f"Per-read logging enabled: {read_log_path}")

    try:
        with pysam.AlignmentFile(in_path, "rb", check_sq=False, threads=threads) as bam_in:

            out_header = update_bam_header(bam_in.header, cmdline=cmdline)

            with pysam.AlignmentFile(out_path, "wb", header=out_header,
                                     threads=threads) as bam_out:

                for rec in bam_in.fetch(until_eof=True):
                    metrics["reads_total"] += 1
                    try:
                        new_rec = remove_barcode_kinetics(rec, drop_mm, metrics, read_logger)
                        bam_out.write(new_rec)
                    except TrimError as exc:
                        metrics["reads_failed"] += 1
                        logging.error("Read %s failed trim: %s", rec.query_name, exc)
                        raise 
    finally:
        # Close the read log file if it was opened
        if read_log_file:
            read_log_file.close()
            logging.info(f"Per-read log written to: {read_log_path}")

    # Summary
    logging.info(
        "Done. Total=%d  fixed=%d  skipped=%d  failed=%d  MM/ML removed=%d",
        metrics["reads_total"],
        metrics["reads_fixed"],
        metrics["reads_skipped"],
        metrics["reads_failed"],
        metrics["mm_ml_removed"],
    )

    if metrics["reads_failed"]:
        raise SystemExit("Some reads failed trimming; see log for details.")

def get_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Trim PacBio kinetics arrays (fi/fp/ri/rp) in cases where lima <2.5 "
                    "was used for demultiplexing. These files retained kinetics values for " 
                    "barcodes but not sequence values which resulted in mismatches."
    )
    p.add_argument("input_bam", help="Input uBAM file (should have been demultiplexed with lima <2.5 and should have kinetics tags")
    p.add_argument("output_bam", help="Output BAM with fixed kinetics")
    p.add_argument("--threads", type=int, default=4, help="htslib threads for reading/writing (default: 4)")
    p.add_argument("--drop_mm", action="store_true", default=False, help="Remove existing MM/ML methylation tags")
    p.add_argument("-v", "--verbose", action="count", default=0, help="Increase log verbosity (-v or -vv)")
    
    return p


def main(argv: list[str] | None = None) -> None:

    args = get_parser().parse_args(argv)

    log_level = logging.WARNING if args.verbose == 0 else \
                logging.INFO if args.verbose == 1 else \
                logging.DEBUG
    
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if log_level == logging.DEBUG: 
        base_name = os.path.splitext(args.output_bam)[0]
        read_log_path = f"{base_name}_read_trim_log.csv"
    else:
        read_log_path = None

    try:
        process_file(args.input_bam, args.output_bam,
                     threads=args.threads,
                     drop_mm=args.drop_mm,
                     read_log_path=read_log_path)
        
    except Exception as exc:
        logging.critical("Failed: %s", exc)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
