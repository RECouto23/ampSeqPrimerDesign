#!/usr/bin/env python
# coding: utf-8

import primer3
import pandas
import pybedtools as pbt
import os
import argparse
import shutil
import yaml
from datetime import datetime

# ── Environment YAML generation ──────────────────────────────────────────────

def write_environment_yaml(out_dir):
    """Write a conda environment YAML for this script's dependencies."""
    env = {
        'name': 'ampseq_primer_design',
        'channels': ['bioconda', 'conda-forge', 'defaults'],
        'dependencies': [
            'python=3.10',
            'primer3-py',
            'pandas',
            'pybedtools',
            'bedtools',
            'pyyaml',
            'pip',
        ]
    }
    yaml_path = os.path.join(out_dir, 'environment.yaml')
    with open(yaml_path, 'w') as fh:
        yaml.dump(env, fh, default_flow_style=False, sort_keys=False)
    print(f'  Environment YAML written → {yaml_path}')

# ── BEDTools PATH resolution ──────────────────────────────────────────────────

def resolve_bedtools():
    """Locate bedtools on PATH; raise a clear error if not found."""
    bt_path = shutil.which('bedtools')
    if bt_path is None:
        raise EnvironmentError(
            'bedtools not found on PATH. Please install bedtools and ensure it '
            'is accessible (e.g. via conda activate or module load).'
        )
    return bt_path

# ── Argument parsing ──────────────────────────────────────────────────────────

now = datetime.today().strftime('%Y%m%d_%H%M%S')

parser = argparse.ArgumentParser(
    prog='Primer3 Design Script',
    description='Submit sites in BED format to generate primer sets using Primer3.'
)
parser.add_argument('-f', '--fileName',
                    help='Input BED file containing sites of targets (tab-delimited, no header).')
parser.add_argument('-o', '--outDir', default='Outputs',
                    help='Output directory for primer CSV files and working files.')
parser.add_argument('-s', '--ampliconSize', default=200,
                    help='Amplicon target size in base pairs.')
parser.add_argument('-b', '--buffer', default=50,
                    help='Acceptable buffer range for amplicon size. '
                         'buffer-ampliconSize < amplicons < buffer+ampliconSize.')
parser.add_argument('-d', '--designSpace', default=1.5,
                    help='Size of design space for primers as a ratio of amplicon size. Must be >1.0.')
parser.add_argument('-t', '--meltingTemperature', default=60,
                    help='Ideal primer melting temperature, in ˚C.')
parser.add_argument('-n', '--numberSets', default=3,
                    help='Number of unique sets of primers to design per locus.')
parser.add_argument('--FWDHandle', default='ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
                    help='Primer handle sequence for FWD primers.')
parser.add_argument('--REVHandle', default='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
                    help='Primer handle sequence for REV primers.')
parser.add_argument('-p', '--probe01', default=0,
                    help='Binary input for whether or not to design a probe. Yes (1), No (0).')
parser.add_argument('-r', '--referencePath',
                    help='Path to reference genome FASTA file used to fetch design-space sequences.')

args = parser.parse_args()

# ── Setup ─────────────────────────────────────────────────────────────────────

print('\nScript to design primers using target sites in BED file format. Relies on Primer3 \n'
      'to do the primer designs. Use the -h parameter for list of input parameters.\n\n')

resolve_bedtools()

isAmpSeq     = True
ampSeqBuffer = int(int(args.ampliconSize) * float(args.designSpace) / 2)
primersOptTM = int(args.meltingTemperature)
numSets      = int(args.numberSets)

# Ask Primer3 for a larger candidate pool so deduplication has room to work.
# Request 5x the desired sets, capped at a reasonable ceiling.
candidatePool = min(numSets * 5, 25)

def resolve_run_dir(base_dir):
    """
    Return a run-specific subdirectory: <base_dir>/YYYYMMDD/RunX
    where X increments based on how many Run* folders already exist for today.
    """
    today    = datetime.today().strftime('%Y%m%d')
    date_dir = os.path.join(base_dir, today)
    os.makedirs(date_dir, exist_ok=True)

    existing = [d for d in os.listdir(date_dir)
                if os.path.isdir(os.path.join(date_dir, d)) and d.startswith('Run')]
    run_num  = len(existing) + 1
    run_dir  = os.path.join(date_dir, f'Run{run_num}')
    os.makedirs(run_dir, exist_ok=True)
    return run_dir

run_dir = resolve_run_dir(str(args.outDir))
print(f'  Output directory -> {run_dir}\n')
write_environment_yaml(run_dir)

inputbed = pandas.read_csv(
    args.fileName, sep='\t', header=None,
    names=['chrom', 'start', 'stop', 'name', 'score', 'strand']
)
inputbed.reset_index(drop=True, inplace=True)

primerSeqs = pandas.DataFrame()

# ── Show input sites ──────────────────────────────────────────────────────────

print('----------DESIGNING PRIMERS FOR THESE SITES----------\n')
print(inputbed)

# ── Expand design space and fetch sequences ───────────────────────────────────

if isAmpSeq:
    for a in range(len(inputbed)):
        inputbed.loc[a, 'start'] = int(inputbed.loc[a, 'start']) - ampSeqBuffer
        inputbed.loc[a, 'stop']  = int(inputbed.loc[a, 'stop'])  + ampSeqBuffer

    inputbed.to_csv(
        os.path.join(run_dir, 'DesignSpaceBED_' + now + '.txt'),
        sep='\t', header=None, index=None
    )

pbtObj = pbt.BedTool(os.path.join(run_dir, 'DesignSpaceBED_' + now + '.txt'))
seq = pbtObj.getfasta(
    fi=args.referencePath,
    bed=os.path.join(run_dir, 'DesignSpaceBED_' + now + '.txt'),
    fo=os.path.join(run_dir, 'DesignSequences_' + now + '.txt')
)

fastaSeqs = pandas.read_csv(
    os.path.join(run_dir, 'DesignSequences_' + now + '.txt'),
    header=None, sep='\t'
)

seqNames = []
dropme   = []
for a in range(0, len(fastaSeqs), 2):
    seqNames.extend(fastaSeqs.loc[a, 0])
    dropme.extend([a])

fastaSeqs = fastaSeqs.drop(dropme).reset_index(drop=True)

# ── Primer design loop ────────────────────────────────────────────────────────

print('\n\n----------BEGINNING PRIMER DESIGN----------\n')

runningList = []

for a in range(0, len(fastaSeqs)):
    print('Designing for: ' + str(inputbed.loc[a, 'name']))

    primers = primer3.bindings.design_primers(
        {
            'SEQUENCE_ID':       inputbed.loc[a, 'name'],
            'SEQUENCE_TEMPLATE': fastaSeqs.values[a][0],
        },
        {
            'PRIMER_PICK_LEFT_PRIMER':        1,
            'PRIMER_PICK_RIGHT_PRIMER':       1,
            'PRIMER_PICK_INTERNAL_OLIGO':     args.probe01,
            'PRIMER_NUM_RETURN':              candidatePool,
            'PRIMER_OPT_SIZE':                20,
            'PRIMER_MIN_SIZE':                17,
            'PRIMER_MAX_SIZE':                25,
            'PRIMER_OPT_TM':                  primersOptTM,
            'PRIMER_MIN_TM':                  primersOptTM - 3,
            'PRIMER_MAX_TM':                  primersOptTM + 3,
            'PRIMER_INTERNAL_OPT_TM':         primersOptTM + 5,
            'PRIMER_INTERNAL_MAX_TM':         (primersOptTM + 5) + 3,
            'PRIMER_INTERNAL_MIN_TM':         (primersOptTM + 5) - 3,
            'PRIMER_MIN_GC':                  20.0,
            'PRIMER_MAX_GC':                  80.0,
            'PRIMER_MAX_POLY_X':              5,
            'PRIMER_SALT_MONOVALENT':         50.0,
            'PRIMER_DNA_CONC':                50.0,
            'PRIMER_MAX_NS_ACCEPTED':         0,
            'PRIMER_MAX_SELF_ANY':            12,
            'PRIMER_MAX_SELF_END':            8,
            'PRIMER_PAIR_MAX_COMPL_ANY':      12,
            'PRIMER_PAIR_MAX_COMPL_END':      8,
            'PRIMER_PRODUCT_SIZE_RANGE':      [
                int(args.ampliconSize) - int(args.buffer),
                int(args.ampliconSize) + int(args.buffer)
            ],
        }
    )

    primerSeqs.at[a, 0] = inputbed.loc[a, 'name']

    # ── Deduplicate: walk the candidate pool, collect unique pairs up to numSets ──
    seen_left  = set()
    seen_right = set()
    uniqueSets = []

    for b in range(0, candidatePool):

        # Stop if Primer3 has no more candidates to offer
        leftReturned  = primers['PRIMER_LEFT_NUM_RETURNED']
        rightReturned = primers['PRIMER_RIGHT_NUM_RETURNED']
        if b >= leftReturned and b >= rightReturned:
            break

        rightSeq = 'No Right Primer Designed'
        rightTm  = 'N/A'
        leftSeq  = 'No Left Primer Designed'
        leftTm   = 'N/A'
        intSeq   = 'No Probe Designed'
        intTm    = 'N/A'
        size     = 'N/A'
        penalty  = 'N/A'

        if b < rightReturned:
            rightSeq = primers['PRIMER_RIGHT_' + str(b) + '_SEQUENCE']
            rightTm  = round(primers['PRIMER_RIGHT_' + str(b) + '_TM'], 1)
        if b < leftReturned:
            leftSeq  = primers['PRIMER_LEFT_'  + str(b) + '_SEQUENCE']
            leftTm   = round(primers['PRIMER_LEFT_'  + str(b) + '_TM'], 1)
        if b < primers['PRIMER_INTERNAL_NUM_RETURNED']:
            intSeq   = primers['PRIMER_INTERNAL_' + str(b) + '_SEQUENCE']
            intTm    = round(primers['PRIMER_INTERNAL_' + str(b) + '_TM'], 1)
        if b < primers.get('PRIMER_PAIR_NUM_RETURNED', 0):
            size    = primers['PRIMER_PAIR_' + str(b) + '_PRODUCT_SIZE']
            penalty = round(primers['PRIMER_PAIR_' + str(b) + '_PENALTY'], 3)

        # Skip if either primer sequence duplicates one already accepted
        leftKey  = leftSeq.upper()  if leftSeq  not in ('No Left Primer Designed',  'N/A') else None
        rightKey = rightSeq.upper() if rightSeq not in ('No Right Primer Designed', 'N/A') else None

        if leftKey in seen_left or rightKey in seen_right:
            continue

        if leftKey:
            seen_left.add(leftKey)
        if rightKey:
            seen_right.add(rightKey)

        uniqueSets.append([leftSeq, leftTm, rightSeq, rightTm, intSeq, intTm, size, penalty])

        if len(uniqueSets) == numSets:
            break

    # Report however many unique sets were found
    for setIdx, s in enumerate(uniqueSets):
        leftSeq, leftTm, rightSeq, rightTm, intSeq, intTm, size, penalty = s
        print('\tSet ' + str(setIdx + 1) + ':'
              + ' FWD - ' + str(leftSeq).upper()
              + ' REV - ' + str(rightSeq).upper())

        runningList.append([
            inputbed.loc[a, 'name'],
            'Set ' + str(setIdx + 1),
            leftSeq, leftTm,
            rightSeq, rightTm,
            intSeq, intTm,
            size, penalty,
        ])

    if len(uniqueSets) < numSets:
        print(f'\t  Warning: only {len(uniqueSets)} unique primer set(s) found for '
              f'{inputbed.loc[a, "name"]} (requested {numSets}).')

print('\n----------PRIMER DESIGN COMPLETE----------')

primerDF = pandas.DataFrame(
    runningList,
    columns=[
        'Site', 'Set Number',
        'Left Primer', 'Left Primer Tm (C)',
        'Right Primer', 'Right Primer Tm (C)',
        'Probe', 'Probe Tm (C)',
        'Product Size', 'Product Penalty Score',
    ]
)

# ── Append handles ────────────────────────────────────────────────────────────

print('\n----------APPENDING HANDLES----------')

primerDF_Adapted = primerDF.copy()

for a in range(len(primerDF_Adapted)):
    primerDF_Adapted.loc[a, 'Left Primer']  = (
        str(args.FWDHandle) + str(primerDF_Adapted.loc[a, 'Left Primer'])
    )
    primerDF_Adapted.loc[a, 'Right Primer'] = (
        str(args.REVHandle) + str(primerDF_Adapted.loc[a, 'Right Primer'])
    )

# ── Write CSVs ────────────────────────────────────────────────────────────────

gsp_path     = os.path.join(run_dir, 'DesignedPrimers_GSP_'     + now + '.csv')
adapted_path = os.path.join(run_dir, 'DesignedPrimers_Adapted_' + now + '.csv')

primerDF.to_csv(gsp_path, index=False)
primerDF_Adapted.to_csv(adapted_path, index=False)

print(f'\n  GSP primers     → {gsp_path}')
print(f'  Adapted primers → {adapted_path}')

print('\n----------COMPLETE----------\n\n\n')
