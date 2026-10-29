"""Prepare public county age/sex cells for adjusted cluster comparisons.
Run from the supplement root; downloads are cached. No raw records are exported.
"""
from pathlib import Path
from collections import Counter
import csv,os,urllib.request,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
PRIMARY=Path(os.environ.get('STARS_SOURCE_DIR') or ROOT.parent)
CACHE=ROOT/'data/external/census';OUT=ROOT/'data/processed/clustering'
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(name,rows):
 with (OUT/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
sources=[];pops={};totals={}
for name,url,yearmap in [
 ('cc-est2020int-alldata-30.csv','https://www2.census.gov/programs-surveys/popest/datasets/2010-2020/intercensal/county/asrh/cc-est2020int-alldata-30.csv',{v+2:v+2010 for v in range(2,10)}),
 ('cc-est2021-alldata-30.csv','https://www2.census.gov/programs-surveys/popest/datasets/2020-2021/counties/asrh/cc-est2021-alldata-30.csv',{2:2020,3:2021})]:
 p=CACHE/name
 if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=60).read())
 sources.append(dict(file=str(p.relative_to(ROOT)),url=url,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 for r in read(p):
  y=yearmap.get(int(r['YEAR']))
  if y is None:continue
  assert r['STATE']=='30' and int(r['SUMLEV'])==50
  fips=r['COUNTY'].zfill(3);age=int(r['AGEGRP'])
  if age==0:totals[fips,y]=int(r['TOT_POP']);continue
  group=0 if age<=3 else (age-2)//2
  # groups: <15, 15-24, 25-34, ..., 75-84, 85+.
  if age==18:group=8
  for sex,col in [('MALE','TOT_MALE'),('FEMALE','TOT_FEMALE')]:
   key=(fips,y,group,sex);pops[key]=pops.get(key,0)+int(r[col])
assert len(totals)==560
context=read(OUT/'county_context.csv');names={r['COUNTY']:r['FIPS'].zfill(3) for r in context};byf={v:k for k,v in names.items()}
for r in read(OUT/'county_year_mortality.csv'):
 key=r['FIPS'].zfill(3),int(r['YEAR'])
 assert totals[key]==int(r['POPULATION']),(key,totals[key],r['POPULATION'])
for (f,y),total in totals.items():assert sum(v for (ff,yy,g,s),v in pops.items() if ff==f and yy==y)==total
stars=read(PRIMARY/'R - Data Exploration/20230831_Montana_Suicide_Database.csv')
def known(x):return x not in ('','NA',None)
def agegroup(age):
 if age<15:return 0
 if age>=85:return 8
 return (age-15)//10+1
agecount=Counter();sexcount=Counter();audit=[]
for f,name in byf.items():
 for y in range(2012,2022):
  d=[r for r in stars if r['DEM_RES03CountyCorrected'].strip()==name and int(r['DEATH_2Year'])==y]
  am=sm=0
  for r in d:
   if not known(r['DEM_AGE1Age']):am+=1;continue
   a=agegroup(int(float(r['DEM_AGE1Age'])));agecount[f,y,a]+=1
   if not known(r['DEM_GEN2GenderNum']):sm+=1;continue
   sex='MALE' if r['DEM_GEN2GenderNum']=='0' else 'FEMALE'
   assert r['DEM_GEN2GenderNum'] in ('0','1')
   sexcount[f,y,a,sex]+=1
  audit.append(dict(FIPS=f,COUNTY=name,YEAR=y,total_deaths=len(d),age_missing=am,
    sex_missing_among_known_age=sm,age_known=len(d)-am,age_sex_known=len(d)-am-sm))
labels=['<15','15-24','25-34','35-44','45-54','55-64','65-74','75-84','85+']
rows=[];sexrows=[]
for f,y in sorted(totals):
 for g,label in enumerate(labels):
  population=sum(pops[f,y,g,s] for s in ['MALE','FEMALE'])
  rows.append(dict(FIPS=f,COUNTY=byf[f],YEAR=y,AGE_GROUP=label,DEATHS=agecount[f,y,g],POPULATION=population))
  if y>=2013:
   for sex in ['MALE','FEMALE']:
    sexrows.append(dict(FIPS=f,COUNTY=byf[f],YEAR=y,AGE_GROUP=label,SEX=sex,
      DEATHS=sexcount[f,y,g,sex],POPULATION=pops[f,y,g,sex]))
assert sum(r['total_deaths'] for r in audit)==2650
assert sum(r['DEATHS'] for r in rows)==sum(r['age_known'] for r in audit)
assert sum(r['DEATHS'] for r in sexrows)==sum(r['age_sex_known'] for r in audit if r['YEAR']>=2013)
assert all(r['POPULATION']>0 or r['DEATHS']==0 for r in rows+sexrows)
write('county_year_age_mortality.csv',rows);write('county_year_age_sex_mortality.csv',sexrows)
write('county_adjusted_coverage.csv',audit);write('county_adjusted_source_manifest.csv',sources)
print('Prepared',len(rows),'age cells and',len(sexrows),'age/sex cells; numerator totals',sum(r['DEATHS'] for r in rows),sum(r['DEATHS'] for r in sexrows))
