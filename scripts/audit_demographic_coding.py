"""Audit the supplied STARS CSV against its dataset key without changing records.

Run from the supplement root:
    python scripts/audit_demographic_coding.py

Uses Python's standard library only. Primary data stay outside the supplement.
"""
import csv
import hashlib
import io
import os
from pathlib import Path
from collections import Counter
import zipfile
import xml.etree.ElementTree as ET

def main():
    root = Path(__file__).resolve().parents[1]
    primary = Path(os.environ.get("STARS_SOURCE_DIR") or root.parent)
    key = primary / "20230528 MT Dataset Key.docx"
    source = primary / "R - Data Exploration/20230831_Montana_Suicide_Database.csv"
    out = root / "outputs/data_quality/coding_audit"
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(key.read_bytes())) as package:
        xml = ET.fromstring(package.read("word/document.xml"))
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    key_rows = []
    for row in xml.findall(".//w:tr",ns):
        cells = ["".join(t.text or "" for t in cell.findall(".//w:t",ns))
                 for cell in row.findall("w:tc",ns)]
        if cells and cells[0].startswith(("DEM_EDU","DEM_MAR","DEM_OCC3","DEM_OCC5")):
            key_rows.append({"field":cells[0],"key_definition":" | ".join(cells[1:])})
    with source.open(encoding="utf-8-sig",newline="") as f:
        rows = list(csv.DictReader(f))
    def missing(x):
        return x is None or str(x).strip() in ("","NA")
    def write(name,records,columns=None):
        if columns is None: columns = list(records[0])
        with (out/name).open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=columns); w.writeheader();w.writerows(records)
    write("key_definitions.csv",key_rows)
    variables = [v for v in rows[0] if v.startswith(("DEM_EDU","DEM_MAR"))]
    yearly = []
    for year in sorted({int(r["DEATH_2Year"]) for r in rows}):
        records=[r for r in rows if int(r["DEATH_2Year"])==year]
        for v in variables:
            counts=Counter("MISSING" if missing(r[v]) else r[v].strip() for r in records)
            known=sum(n for label,n in counts.items() if label!="MISSING")
            for label,n in sorted(counts.items()):
                yearly.append(dict(year=year,field=v,category=label,n=n,year_records=len(records),
                    known=known,percent_of_year=100*n/len(records),
                    percent_of_known=100*n/known if known and label!="MISSING" else ""))
    write("categories_by_year.csv",yearly)
    mapping={1:"NO HIGH SCHOOL",2:"SOME HIGH SCHOOL",3:"HIGH SCHOOL OR GED",
             4:"SOME COLLEGE",5:"ASSOCIATES",6:"BACHELORS",7:"MASTERS",8:"DOCTORATE"}
    mismatch=[]; pairs=Counter()
    marital_pairs=Counter()
    for row in rows:
        year=int(row["DEATH_2Year"])
        text=row["DEM_EDU1EducationAchieved"].strip()
        number=row["DEM_EDU1EducationAchievedNum"].strip()
        pairs[(year,text or "MISSING",number or "MISSING")]+=1
        if not missing(text) and not missing(number):
            if mapping.get(int(float(number)))!=text:
                mismatch.append({"year":year,"text":text,"number":number})
        elif missing(text)!=missing(number):
            mismatch.append({"year":year,"text":text,"number":number})
        marital_pairs[(year,row["DEM_MAR1MaritalStatus"] or "MISSING",
            row["DEM_MAR2EverMarried"],row["DEM_MAR3Married"])]+=1
    write("education_text_numeric_pairs.csv",
          [dict(year=y,text=t,numeric_code=n,count=c) for (y,t,n),c in sorted(pairs.items())])
    write("education_mismatches.csv",mismatch,["year","text","number"])
    write("marital_redundant_field_pairs.csv",
          [dict(year=y,status=t,ever_married=e,married=m,count=c)
           for (y,t,e,m),c in sorted(marital_pairs.items())])
    binary=[]
    for v in ["DEM_OCC3RancherFarmer","DEM_OCC5ArmedForcesEver"]:
        for value,n in sorted(Counter(r[v] for r in rows).items()):
            binary.append(dict(field=v,value=value,count=n))
    write("observed_binary_codes.csv",binary)
    summaries=[]
    for year in range(2012,2022):
        z=[r for r in rows if int(r["DEATH_2Year"])==year]
        educ=Counter(r["DEM_EDU1EducationAchieved"] for r in z)
        marital=Counter(r["DEM_MAR1MaritalStatus"] for r in z)
        summaries.append(dict(year=year,records=len(z),
            education_known=sum(n for k,n in educ.items() if not missing(k)),
            some_high_school=educ["SOME HIGH SCHOOL"],high_school_ged=educ["HIGH SCHOOL OR GED"],
            college_degree=sum(educ[k] for k in ["ASSOCIATES","BACHELORS","MASTERS","DOCTORATE"]),
            marital_known=sum(n for k,n in marital.items() if not missing(k)),
            divorced=marital["DIVORCED"],never_married=marital["NEVER MARRIED"],single=marital["SINGLE"],
            married_not_separated=sum(n for k,n in marital.items()
                if k.startswith("MARRIED") and "NOT SEPARATED" in k),
            married_separated=marital["MARRIED (SEPARATED)"],
            widowed=marital["WIDOWED"]))
    write("annual_coding_summary.csv",summaries)
    write("input_checksums.csv",[dict(input=str(f.relative_to(primary)),
        sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in [key,source]])
    checks=[dict(check="2764 input records",passed=len(rows)==2764),
        dict(check="all ten study years",passed={int(r["DEATH_2Year"]) for r in rows}==set(range(2012,2022))),
        dict(check="education label numeric consistency",passed=not mismatch),
        dict(check="education frequency accounting",passed=sum(c for _,c in pairs.items())==len(rows))]
    write("audit_checks.csv",checks)
    report="""Dataset key and CSV demographic coding audit

Sources: 20230528 MT Dataset Key.docx and
R - Data Exploration/20230831_Montana_Suicide_Database.csv.
This audit does not change primary records or regression recodes.

Findings
- The key explicitly identifies marital status and ever-married as UNCLEAR.
- It does not document year-specific coding crosswalks or identify which annual
  categories should be merged, relabeled or treated as unknown.
- Education numeric codes and text labels agree for all available paired values.
  They encode the same information; neither supplies independent verification.
- There are no HIGH SCHOOL OR GED records in 2013-2017 or 2019 despite known
  education in those years. The distinction from SOME HIGH SCHOOL cannot be
  reconstructed reliably from these two fields.
- No married records are labeled in 2013-2016. In 2020 and 2021, every record
  labeled as married is labeled MARRIED (SEPARATED), with no not-separated category.
  These patterns strongly question comparability but do not prove the intended
  replacement category.
- Redundant ever-married and married fields reflect the same status categories.
  They do not establish whether separated labels should instead mean married.
- The key's farm/veteran coding describes 1/2, whereas the supplied CSV uses 0/1.
  The regression implementation follows observed 0/1 coding. The key cannot
  be treated as an exact specification of every export field.

Interpretation
The degree indicator combines all below-associate attainment labels. Partial
attainment retains the recorded SOME HIGH SCHOOL and SOME COLLEGE distinction.
Both classifications are analyzed; neither validates the underlying labels.
Divorced and separated are combined in the principal regressions. A same-sample
comparison omitting marital status assesses its influence on other estimates.
Year effects account for calendar differences in the intercept but do not
correct exposure misclassification. Frequency patterns alone do not identify
replacement labels or establish the reason for absent categories.

The audit establishes internal consistency and recording discontinuities. It
cannot determine the true marital status or educational history of individuals.
"""
    (out/"audit_report.txt").write_text(report,encoding="utf-8")
    print(report)
    print("Education mismatches:",len(mismatch))
    print("Results:",out)
    if not all(c["passed"] for c in checks):raise SystemExit("Audit checks failed")

if __name__=="__main__":
    main()
