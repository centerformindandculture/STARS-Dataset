"""Build county comparison county profiles and a FIPS-joined cluster map.

Run after Rscript scripts/run_clustering.R. Uses bundled public aggregate data
and cached 2020 Census boundaries; no private records or network needed.
Python standard library suffices for CSV, HTML and SVG. --png needs Pillow.
"""
from pathlib import Path
import os, argparse, csv, html, json, math, statistics, sys
from collections import defaultdict

COLORS = {1:'#E69F00', 2:'#56B4E9', 3:'#009E73'}
LABELS = [
 ('MEDAGE','Median age, years'),
 ('AIAN_ALONE_PCT','American Indian / Alaska Native alone, %'),
 ('MEDINCOME','Median household income, 2011 dollars'),
 ('POVERTY_PCT','Population below poverty threshold, %'),
 ('SINGLE_HH_PCT','One-person households, %'),
 ('DIVSEP_PCT','Divorced or separated, age 15+, %'),
 ('DEGREE_PCT','Associate degree or higher, age 25+, %'),
 ('UNINSURED_PCT','Uninsured population, %')]
NOTES = [
 'Contextual characteristics are unweighted means (SD) across counties. Each county contributes equally; these are not pooled resident percentages. Seven characteristics use 2007-2011 ACS five-year estimates; uninsured percentage uses 2008-2012 ACS.',
 'Annual deaths per county are summarized using each county\'s 2012-2021 mean annual count; SD is across counties. Total deaths and person-years describe aggregate burden and exposure.',
 'Mortality rates pool deaths and annual populations over 2012-2021, per 100,000 person-years. Parentheses are descriptive exact Poisson 95% intervals with fixed population denominators. Rates are crude and are not age-standardized.',
 'P-values are omnibus comparisons across the three clusters from Poisson models with categorical year and county-clustered sandwich with finite-sample correction covariance (F reference, 2 and 55 df). The rate model includes log(population) as an offset; the count model does not. Displayed p-values are Holm-adjusted across the four preferred county-count, county-rate, month-rate and weekday-rate tests. Raw and adjusted values are retained in the numeric export.',
 'No p-values are given for defining cluster characteristics, because clustering was designed to separate counties on these measures. Totals are descriptive; their differing number of counties is not itself tested. Census uncertainty and cluster-selection uncertainty are not propagated.',
 'Adjusted rates use Poisson models with cluster, categorical year and age or joint age-by-sex strata, with population offsets. They are standardized to a common pooled Montana demographic/year distribution. Age bands are <15, 15-24, 25-34, 35-44, 45-54, 55-64, 65-74, 75-84 and 85+. Intervals use county-clustered sandwich with finite-sample correction covariance and t(55) limits.',
 'Age adjustment covers 2012-2021 (2,650 county-matched deaths). Age/sex adjustment covers 2013-2021 (2,441 deaths), because sex is unavailable in 2012. These are internally standardized rates; the two periods differ. Matched-period crude rates are retained in county_k3_adjusted_rates.csv.',
 'The two demographic-standardization tests form a separate Holm-adjusted family. Population-estimate and cluster-selection uncertainty are not propagated.',
 'All 56 counties are included, including Petroleum with zero deaths. The analysis covers 2,650 county-matched deaths; 114 records without a usable Montana county are excluded.']

def read(p):
 with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def write(p, rows):
 with p.open('w',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def fmt_p(p):return '<0.0001' if p<.0001 else f'{p:.4f}'
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--outputs',default=str(Path(os.environ.get('STARS_OUTPUT_DIR','outputs'))/'clustering/geographic'));ap.add_argument('--png',action='store_true');a=ap.parse_args()
 root=Path(__file__).resolve().parents[1];out=root/a.outputs
 context=read(root/'data/processed/clustering/county_context.csv')
 assignment=read(out/'county_k3_assignments.csv'); lookup={r['COUNTY']:int(r['CLUSTER']) for r in assignment}
 assert len(lookup)==56 and set(lookup)=={r['COUNTY'] for r in context}
 assert [list(lookup.values()).count(g) for g in (1,2,3)]==[28,6,22]
 mortality={int(r['cluster']):r for r in read(out/'county_k3_mortality.csv')}
 tests={r['metric']:r for r in read(out.parent/'overview/mortality_omnibus_tests.csv') if r['domain']=='county' and r['k']=='3'}
 year=read(out/'county_k3_county_year_analysis.csv'); bycounty=defaultdict(list)
 for r in year:bycounty[r['COUNTY']].append(r)
 assert len(year)==560 and all(len(v)==10 for v in bycounty.values())
 assert sum(int(r['DEATHS']) for r in year)==2650
 assert all(sum(int(x['DEATHS']) for x in bycounty[r['COUNTY']])==0 for r in context if r['COUNTY']=='Petroleum')
 rows=[];numeric=[]
 def add(label,values,p='',raw='',stat='descriptive',means=None,sds=None,p_exact=''):
  rows.append(dict(Characteristic=label,**{f'Cluster {g}':values[g-1] for g in (1,2,3)},**{'Mortality p (Holm)':p}))
  for g in (1,2,3):numeric.append(dict(characteristic=label,cluster=g,display=values[g-1],statistic=stat,estimate='' if means is None else means[g-1],sd='' if sds is None else sds[g-1],p_raw=raw,p_holm=p_exact))
 add('Counties, n',[str(list(lookup.values()).count(g)) for g in (1,2,3)])
 for key,label in LABELS:
  groups=[[float(r[key]) for r in context if lookup[r['COUNTY']]==g] for g in (1,2,3)]
  means=[statistics.mean(v) for v in groups];sds=[statistics.stdev(v) for v in groups]
  if key=='MEDINCOME':values=[f'{m:,.0f} ({s:,.0f})' for m,s in zip(means,sds)]
  else:values=[f'{m:.1f} ({s:.1f})' for m,s in zip(means,sds)]
  add(label,values,stat='unweighted county mean (SD)',means=means,sds=sds)
 add('Person-years, 2012-2021',[f"{float(mortality[g]['person_years']):,.0f}" for g in (1,2,3)])
 add('Suicide deaths, total',[f"{int(mortality[g]['DEATHS']):,}" for g in (1,2,3)])
 groups=[[sum(int(r['DEATHS']) for r in v)/10 for c,v in bycounty.items() if lookup[c]==g] for g in (1,2,3)]
 means=[statistics.mean(v) for v in groups];sds=[statistics.stdev(v) for v in groups]
 for g in (1,2,3):assert math.isclose(means[g-1],float(mortality[g]['mean_annual_deaths_per_county']),rel_tol=1e-12)
 add('Annual suicide deaths per county, mean (SD)',[f'{m:.2f} ({s:.2f})' for m,s in zip(means,sds)],fmt_p(float(tests['count']['p_holm_preferred'])),tests['count']['p'],'county mean annual count (SD)',means,sds,p_exact=tests['count']['p_holm_preferred'])
 values=[];rates=[]
 for g in (1,2,3):
  r=mortality[g];rate=int(r['DEATHS'])/float(r['person_years'])*1e5;rates.append(rate)
  assert math.isclose(rate,float(r['rate_per_100000']),rel_tol=1e-12)
  values.append(f"{rate:.2f} ({float(r['rate_ci_lower']):.2f}-{float(r['rate_ci_upper']):.2f})")
 add('Suicide deaths / 100,000 person-years (95% CI)',values,fmt_p(float(tests['rate']['p_holm_preferred'])),tests['rate']['p'],'pooled crude rate (95% Poisson CI)',rates,p_exact=tests['rate']['p_holm_preferred'])
 adjusted=read(out/'county_k3_adjusted_rates.csv')
 adjusted_tests={r['adjustment']:r for r in read(out/'county_k3_adjusted_tests.csv')}
 for kind,label in [('age','Age-adjusted deaths / 100,000, 2012-2021 (95% CI)'),('age_sex','Age/sex-adjusted deaths / 100,000, 2013-2021 (95% CI)')]:
  rr=sorted([r for r in adjusted if r['adjustment']==kind],key=lambda r:int(r['cluster']))
  assert [int(r['cluster']) for r in rr]==[1,2,3]
  vals=[f"{float(r['standardized_rate_per_100000']):.2f} ({float(r['ci_lower']):.2f}-{float(r['ci_upper']):.2f})" for r in rr]
  test=adjusted_tests[kind]
  add(label,vals,fmt_p(float(test['p_holm'])),test['p'],
      'model-standardized rate (95% county-clustered CI)',[float(r['standardized_rate_per_100000']) for r in rr],p_exact=test['p_holm'])
 write(out/'county_k3_comparison_table.csv',rows);write(out/'county_k3_comparison_numeric.csv',numeric)
 # Attach FIPS from the same public county context, not county-name geometry matching.
 fips={str(r['FIPS']).zfill(3):dict(county=r['COUNTY'],cluster=lookup[r['COUNTY']]) for r in context}
 boundary_path=root/'data/external/census/montana_counties_2020_5070.json'
 boundaries=json.loads(boundary_path.read_text(encoding='utf-8'));features=boundaries['features']
 assert boundaries['spatialReference']['wkid']==5070
 assert len(features)==56 and {f['attributes']['GEOID'] for f in features}=={'30'+v for v in fips}
 join=[]
 for f in features:
  att=f['attributes'];item=fips[att['COUNTY']]
  assert att['GEOID']=='30'+att['COUNTY']
  # A second name check catches an accidental FIPS lookup mismatch.
  assert att['NAME'].removesuffix(' County').replace(' ','').lower()==item['county'].replace(' ','').lower()
  f['cluster']=item['cluster'];f['county']=item['county']
  join.append(dict(GEOID=att['GEOID'],COUNTY=item['county'],CLUSTER=item['cluster'],COLOR=COLORS[item['cluster']]))
 write(out/'county_k3_map_join.csv',sorted(join,key=lambda r:r['GEOID']))
 points=[p for f in features for ring in f['geometry']['rings'] for p in ring]
 xmin=min(p[0] for p in points);xmax=max(p[0] for p in points);ymin=min(p[1] for p in points);ymax=max(p[1] for p in points)
 W,H=1800,1160;scale=min(1680/(xmax-xmin),850/(ymax-ymin));left=(W-(xmax-xmin)*scale)/2;top=150
 def xy(p):return left+(p[0]-xmin)*scale,top+(ymax-p[1])*scale
 svg=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-label="Montana counties colored by three contextual clusters">', '<rect width="100%" height="100%" fill="white"/>', '<g font-family="Arial, sans-serif" fill="#172B36">', '<text x="60" y="65" font-size="38" font-weight="bold">Montana county clusters</text>', '<text x="60" y="108" font-size="23">Three clusters based on county context; all 56 counties included</text>']
 for f in features:
  d=' '.join('M '+' L '.join(f'{x:.2f},{y:.2f}' for x,y in map(xy,ring))+' Z' for ring in f['geometry']['rings'])
  svg.append(f'<path d="{d}" fill="{COLORS[f["cluster"]]}" fill-rule="evenodd" stroke="white" stroke-width="1.8"><title>{html.escape(f["county"])} County - Cluster {f["cluster"]}</title></path>')
 for g in (1,2,3):
  x=250+(g-1)*480;svg.extend([f'<rect x="{x}" y="1014" width="34" height="28" fill="{COLORS[g]}"/>',f'<text x="{x+48}" y="1038" font-size="25">Cluster {g} - {list(lookup.values()).count(g)} counties</text>'])
 svg.extend(['<text x="60" y="1092" font-size="19">Boundary source: U.S. Census Bureau, 2020 TIGERweb; NAD83 / Conus Albers (EPSG:5070).</text>','<text x="60" y="1124" font-size="19">Colors identify contextual groups and do not represent suicide mortality levels.</text>','</g></svg>'])
 svgtext='\n'.join(svg);(out/'county_k3_map.svg').write_text(svgtext,encoding='utf-8')
 headers=list(rows[0]);table='<table><thead><tr>'+''.join(f'<th>{html.escape(c)}</th>' for c in headers)+'</tr></thead><tbody>'
 for i,r in enumerate(rows):
  table+='<tr'+(' class="mortality"' if i>=11 else '')+'>'+''.join(f'<{ "th" if j==0 else "td" }>{html.escape(str(r[c]))}</{ "th" if j==0 else "td" }>' for j,c in enumerate(headers))+'</tr>'
 table+='</tbody></table>'
 document='<!doctype html><html lang="en"><meta charset="utf-8"><title>Montana county cluster profiles</title><style>body{font:16px Arial,sans-serif;color:#172b36;max-width:1280px;margin:45px auto;padding:0 28px}h1{font-size:28px}table{border-collapse:collapse;width:100%;font-size:15px}th,td{padding:13px 12px;border-bottom:1px solid #dce3e8;text-align:right}th:first-child{text-align:left;width:45%}thead{background:#eef3f5}thead th:nth-child(2){border-top:6px solid #e69f00}thead th:nth-child(3){border-top:6px solid #56b4e9}thead th:nth-child(4){border-top:6px solid #009e73}tbody th{font-weight:400}.mortality{background:#f3f7f8}.notes{font-size:13px;line-height:1.5}svg{width:100%;height:auto}@media print{body{margin:0}table{font-size:10pt}svg{break-before:page}}</style><h1>County cluster comparison profiles</h1><p>Contextual means (SD), with subsequent suicide mortality comparisons, 2012-2021.</p>'+table+'<div class="notes">'+''.join('<p>'+html.escape(n)+'</p>' for n in NOTES)+'</div>'+svgtext+'</html>'
 (out/'county_k3_profiles.html').write_text(document,encoding='utf-8')
 (out/'county_k3_comparison_notes.txt').write_text('\n\n'.join(NOTES),encoding='utf-8')
 if a.png:
  from PIL import Image,ImageDraw,ImageFont
  im=Image.new('RGB',(W,H),'white');dr=ImageDraw.Draw(im)
  def font(n,bold=False):
   for name in (('arialbd.ttf','DejaVuSans-Bold.ttf') if bold else ('arial.ttf','DejaVuSans.ttf')):
    try:return ImageFont.truetype(name,n)
    except OSError:pass
   return ImageFont.load_default(size=n)
  dr.text((60,28),'Montana county clusters',fill='#172B36',font=font(38,True));dr.text((60,80),'Three clusters based on county context; all 56 counties included',fill='#172B36',font=font(23))
  for f in features:
   rings=sorted(f['geometry']['rings'],key=lambda r:abs(sum(r[i][0]*r[i+1][1]-r[i+1][0]*r[i][1] for i in range(len(r)-1))),reverse=True)
   mask=Image.new('L',(W,H));md=ImageDraw.Draw(mask)
   for ring in rings:
    area=sum(ring[i][0]*ring[i+1][1]-ring[i+1][0]*ring[i][1] for i in range(len(ring)-1))
    md.polygon([xy(p) for p in ring],fill=255 if area<0 else 0)
   im.paste(COLORS[f['cluster']],mask=mask)
   for ring in rings:dr.line([xy(p) for p in ring],fill='white',width=2)
  for g in (1,2,3):
   x=250+(g-1)*480;dr.rectangle((x,1014,x+34,1042),fill=COLORS[g]);dr.text((x+48,1014),f'Cluster {g} - {list(lookup.values()).count(g)} counties',fill='#172B36',font=font(25))
  dr.text((60,1074),'Boundary source: U.S. Census Bureau, 2020 TIGERweb; NAD83 / Conus Albers (EPSG:5070).',fill='#172B36',font=font(19));dr.text((60,1106),'Colors identify contextual groups and do not represent suicide mortality levels.',fill='#172B36',font=font(19))
  im.save(out/'county_k3_map.png',dpi=(300,300))
 print('Created comparison table, numeric export, HTML report, SVG map and FIPS join in',out)
 print('Validated 56 FIPS joins, 560 county-years, 2,650 deaths, population-rate accounting and count summaries.')
if __name__=='__main__':main()
