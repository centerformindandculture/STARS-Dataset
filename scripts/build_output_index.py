"""Create linked and CSV inventories of generated analysis outputs."""
from pathlib import Path
from urllib.parse import quote
import csv, os

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/Path(os.environ.get('STARS_OUTPUT_DIR','outputs'))
DESCRIPTIONS={
 'descriptive':'Source-cohort characteristics, recorded categories and missingness.',
 'regressions/principal':'Four sequential adult resident models, with and without race/ethnicity.',
 'regressions/partial_education':'Four-model sequences using partial high-school/college attainment.',
 'regressions/education_comparison':'Joint education indicators, incremental Wald tests and three-group contrasts.',
 'regressions/age_reference_45_54':'Principal and partial-attainment sequences relative to ages 45–54.',
 'regressions/diagnostics':'Input coding, eligibility, age diagnostics, fit checks and covariance-method comparisons.',
 'regressions/sensitivity/age_sex_interaction':'Joint age-by-sex tests, within-sex age contrasts and cell counts.',
 'regressions/sensitivity/marital_adjustment':'Models omitting and including marital status on identical records.',
 'regressions/sensitivity/marital_three_categories':'Four-model sequences separating divorced and separated records, with equality contrasts.',
 'regressions/sensitivity/late_years':'Models restricted to 2017–2019 and 2021 without requiring race.',
 'regressions/sensitivity/age_25_plus':'Models restricted to adults aged 25+, reference age 25–34.',
 'regressions/sensitivity/exclude_2021':'Fully adjusted principal model using supported years 2013–2019.',
 'regressions/sensitivity/sample_selection':'Retained/excluded profiles and minimal models on different eligible samples.',
 'regressions/sensitivity/race_adjustment':'Model without race on the race-complete sample.',
 'clustering/geographic':'County assignments, profiles, maps, counts and crude/standardized mortality rates.',
 'clustering/monthly':'Month assignments, clustering diagnostics and mortality comparisons.',
 'clustering/weekday':'Weekday assignments, clustering diagnostics and mortality comparisons.',
 'clustering/overview':'Joint cluster selection, stability, model diagnostics and multiplicity summaries.',
 'data_quality/comparison_2020':'2020 distributions and missingness compared with other periods.',
 'data_quality/coding_audit':'Dataset-key definitions and annual category recording.',
 'validation':'Cohort, specification and numerical consistency checks.'}

def describe(p):
 if '_replicate_coefficients' in p.stem:return 'Annual ACS replicate coefficient vectors'
 if '_covariance' in p.stem:return 'Coefficient covariance matrix'
 if p.stem.endswith('_long'):return 'Numeric estimates, confidence limits and exact p-values'
 if p.stem.endswith('_notes'):return 'Definitions, reference groups and inference notes'
 if p.stem.endswith('_summary'):return 'Fitted model summary'
 return p.stem.replace('_',' ').capitalize()

def main():
 files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name not in {'INDEX.md','file_inventory.csv'}
              and not set(p.relative_to(OUT).parts)&{'local','reproduced'})
 groups={};rows=[]
 for p in files:
  rel=p.relative_to(OUT).as_posix()
  candidates=[k for k in DESCRIPTIONS if rel.startswith(k+'/')]
  group=max(candidates,key=len) if candidates else 'run_metadata'
  groups.setdefault(group,[]).append(p)
  rows.append(dict(analysis=group,path=rel,description=describe(p),format=p.suffix.lstrip('.'),bytes=p.stat().st_size))
 with (OUT/'file_inventory.csv').open('w',encoding='utf-8',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=['analysis','path','description','format','bytes']);writer.writeheader();writer.writerows(rows)
 lines=['# Output index','','See the [analysis index](../docs/ANALYSIS_INDEX.md) for methods and primary result files.',
        'The [CSV inventory](file_inventory.csv) lists every generated result.','']
 for group,members in sorted(groups.items()):
  lines+=['## '+group.replace('_',' '),'',DESCRIPTIONS.get(group,'Run-level metadata.'),'']
  for details,label in [(False,'Result files'),(True,'Replicate coefficients and covariance matrices')]:
   selected=[p for p in members if ('replicate_details' in p.parts)==details]
   if not selected:continue
   lines+=['<details>','<summary>'+label+'</summary>','','| File | Contents |','| --- | --- |']
   for p in selected:lines.append(f'| [{p.name}]({quote(p.relative_to(OUT).as_posix())}) | {describe(p)} |')
   lines+=['','</details>','']
 (OUT/'INDEX.md').write_text('\n'.join(lines),encoding='utf-8')
 print(f'Indexed {len(files)} files in {len(groups)} analysis groups.')

if __name__=='__main__':main()
