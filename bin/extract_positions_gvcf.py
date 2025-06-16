#!/usr/bin/env python

import sys
import argparse
import logging
import array as arr
from pathlib import Path
from pysam import FastaFile, VariantFile

logger = logging.getLogger()


def read_pos(infile: Path):
    positions = []
    with open(infile, 'r') as inhandle:
        for line in inhandle:
            fields = line.strip().split(':')
            if len(fields) != 2:
                raise Exception(f"Invalid format of genomic position {line}!")
            positions.append((str(fields[0]), int(fields[1])))
    return positions


def process_gvcf(
        vcf_file: Path,
        fasta_file: Path,
        positions: list,
        outfile: Path
        ):
    with VariantFile(str(vcf_file)) as vcf_in, \
        FastaFile(str(fasta_file)) as fasta_in, \
        VariantFile(str(outfile), mode='w', header=vcf_in.header) as vcf_out:
        processed = 0
        for chr, pos in positions:
            refbase = fasta_in.fetch(reference=chr, start=pos-1, end=pos)
            recs = list(vcf_in.fetch(contig=chr, start=pos-1, stop=pos))
            if len(recs) != 1:
                raise Exception(f"Position {chr}:{pos} not found in gVCF file!")
            rec = recs[0]
            if not rec.pos == pos:
                rec.pos = pos
                rec.alleles = (refbase, rec.alleles[1])
            vcf_out.write(rec)
            processed += 1
        logger.info(f"Processed {processed} positions of gVCF file {vcf_file}.")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract positions of interest from gVCF file.",
        epilog="Example: python extract_positions_gvcf.py input.gvcf reference.fa positions.txt poi.gvcf",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "gvcf",
        type=Path,
        help="gVCF file."
    )
    parser.add_argument(
        "fasta",
        type=Path,
        help="Reference fasta file."
    )
    parser.add_argument(
        "positions",
        type=Path,
        help="Text files with genomic positions."
    )
    parser.add_argument(
        "outfile",
        type=Path,
        help="Output vcf file."
    )
    parser.add_argument(
        "-i",
        "--interval",
        type=int,
        help="The desired interval of progress updates.",
        default=1000
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING"
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    positions = read_pos(args.positions)
    process_gvcf(args.gvcf, args.fasta, positions, args.outfile)


if __name__ == "__main__":
    sys.exit(main())
