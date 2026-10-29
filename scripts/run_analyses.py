"""Run analyses, county figures, validation and the output index."""
from pathlib import Path
import argparse, os, shutil, subprocess, sys

ROOT=Path(__file__).resolve().parents[1]

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('mode',choices=['aggregate','full'],nargs='?',default='aggregate')
 parser.add_argument('--rscript',default=shutil.which('Rscript'),help='Path to Rscript executable')
 parser.add_argument('--prepare',action='store_true',help='Rebuild clustering inputs from primary records and workbooks')
 args=parser.parse_args()
 if not args.rscript:parser.error('Rscript is not on PATH; supply --rscript')
 commands=[]
 if args.prepare:
  commands += [[sys.executable,'scripts/prepare_clustering_inputs.py'],
               [sys.executable,'scripts/prepare_county_adjusted.py']]
 commands += [[args.rscript,'scripts/reproduce.R',args.mode]]
 if args.mode=='full':
  commands += [[sys.executable,'scripts/audit_demographic_coding.py'],
               [args.rscript,'scripts/export_descriptive_aggregates.R']]
 commands += [[sys.executable,'scripts/build_county_profiles.py','--png'],
              [sys.executable,'scripts/build_county_map_legend_only.py'],
              [args.rscript,'scripts/validate.R',args.mode],
              [sys.executable,'scripts/build_output_index.py']]
 for command in commands:
  print('Running:',subprocess.list2cmdline(command),flush=True)
  subprocess.run(command,cwd=ROOT,check=True)
 print('Results:',os.environ.get('STARS_OUTPUT_DIR','outputs'))

if __name__=='__main__':main()
