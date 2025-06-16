#!/usr/bin/env python

import sys
import argparse
import logging
import csv
import gzip
import array as arr
from pathlib import Path
from pysam import VariantFile

logger = logging.getLogger()


class VariantSite:
    def __init__(
            self,
            chrom: str,
            pos: int,
            qual: float,
            samples: list[str],
            depths: list[int] | None,
            alleles: list[int] | None,
            phased: list[int] | None):
        self.chrom = chrom
        self.pos = pos
        self.qual = qual
        self.depths = arr.array('L', depths) if depths is not None else arr.array('H', [0] * len(samples))
        self.alleles = arr.array('b', alleles) if alleles is not None else arr.array('b', [-1] * 2 * len(samples))
        self.phased = arr.array('b', phased) if phased is not None else arr.array('b', [False] * len(samples))


def read_vcf(
        vcf_file: Path,
        samples: list[str] | None,
        interval: int=1000
        ):
    variants = {}
    with VariantFile(str(vcf_file)) as vcf_in:
        samples = list(vcf_in.header.samples) if samples is None else samples
        processed = 0
        for rec in vcf_in:
            depths = [rec.samples[sample]["DP"] if rec.samples[sample]["DP"] is not None else 0 for sample in samples]
            alleles = [al if al is not None else -1 for sample in samples for al in rec.samples[sample]["GT"]]
            phased = [rec.samples[sample].phased for sample in samples]
#            variants[(rec.chrom, rec.pos)] = VariantSite(rec.chrom, rec.pos, rec.qual, samples, depths, alleles, phased)
            for eff in rec.id.split(';'):
                variants[eff] = VariantSite(rec.chrom, rec.pos, rec.qual, samples, depths, alleles, phased)
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of VCF file {vcf_file}.")
    return samples, variants


def read_variant_table(infile: Path, reference: str | None=None, interval: int=1000):
    variants = {}
    with open(infile, 'r') as inhandle:
        processed = 0
        for line in inhandle:
            fields = line.strip().split("\t")
            vid = f"{reference}:{fields[0]}" if reference else fields[0]
            variants[vid] = fields[6]
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of input file {infile}.")
    return variants


def add_genotypes_to_annotation(
        annot_file: Path,
        variants: dict,
        samples: list[str],
        output_file: Path,
        reference: str,
        db: dict[str],
        mindp: int=1,
        interval: int=1000
        ):
    open_fn = gzip.open if str(annot_file).endswith('.gz') else open
    with open_fn(annot_file, 'rt') as inhandle, open(output_file, 'w') as outhandle:
        processed = 0
        for line in inhandle:
            line = line.strip()
            if line.startswith('##'):
                print(line, file=outhandle)
                continue
            fields = line.split("\t")            
            if line.startswith('#'):
                print(line, "db", "\t".join(samples), sep="\t", file=outhandle)
                cols = { col: i for i, col in enumerate(fields) }
                continue
            chrom = fields[1].split(":")[0]
            pos = int(fields[1].split(":")[1].split("-")[0])
            if fields[cols['Feature']] != reference: continue
            db_var = '-'
            if db:
                try:
                    vid = fields[cols['HGVSc']].strip()
                    if vid != '-':
                        db_var = db[vid]
                except:
                    pass
                else:
                    if vid != '-':
                        logger.info(f"Found variant {vid} in database with value {db_var}.")
            try:
                rec = variants[fields[0]]
#                rec = variants[(chrom, pos)]
            except:
                logger.warning(f"Couldn't find variant {fields[0]} at position {chrom}:{pos} in VCF file!")
            if rec:
                alleles = (str(al) if al >= 0 else '.' for al in rec.alleles)
                gts = [ (al1, al2) if dp >= mindp else ('.', '.') for al1, al2, dp in zip(alleles, alleles, rec.depths) ]
                gts_phased = [ f"{al1}|{al2}" if phased else f"{al1}/{al2}" for (al1, al2), phased in zip(gts, rec.phased) ]
            else:
                gts_phased = ['./.'] * len(samples)
            print(line, db_var, "\t".join(gts_phased), sep="\t", file=outhandle)
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of annotation file {annot_file}.")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Add genotypes to VEP annotation report.",
        epilog="Example: python add_genotypes_to_annotation.py annot.tab input.vcf annot_with_genotypes.tab",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "annot",
        type=Path,
        help="VEP annotation report in tab format."
    )
    parser.add_argument(
        "vcf",
        type=Path,
        help="Multisample VCF file with phased genotypes."
    )
    parser.add_argument(
        "outfile",
        type=Path,
        help="Output file for annotation report with genotypes in tab format."
    )
    parser.add_argument(
        "--variants",
        type=Path,
        required=False,
        help="Additional variant annotation."
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        required=False,
        help="Reference id to use for annotations."
    )
    parser.add_argument(
        "-s",
        "--samples",
        nargs='+',
        type=str,
        required=False,
        help="Names of samples to add genotypes."
    )
    parser.add_argument(
        "-d",
        "--mindepth",
        type=int,
        help="Minium depth needed to be considered as valid genotype call.",
        default=1
    )
    parser.add_argument(
        "--only_biallelic",
        help="Use only biallelic variants.",
        action='store_true'
    )
    parser.add_argument(
        "--include_indels",
        help="Include INDEL variants.",
        action='store_true'
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

    samples, variants = read_vcf(args.vcf, args.samples, args.interval)
#    reference = "NM_000492.4"
    if args.variants:
        db = read_variant_table(args.variants, args.reference, args.interval)
    else:
        db = None

    add_genotypes_to_annotation(args.annot, variants, samples, args.outfile, args.reference, db, args.mindepth, args.interval)

if __name__ == "__main__":
    sys.exit(main())
