"""Export month and weekday context from the study's external workbooks."""
from pathlib import Path
import os, csv, hashlib, calendar
from openpyxl import load_workbook
ROOT=Path(__file__).resolve().parents[1]
PRIMARY=Path(os.environ.get("STARS_SOURCE_DIR") or ROOT.parent)
OUT=ROOT/"data/processed/clustering"
OUT.mkdir(parents=True,exist_ok=True)
sources=[]
def read(relative, columns):
    path=PRIMARY/relative
    book=load_workbook(path,read_only=True,data_only=True)
    try:
        rows=iter(book["Data"].values);header=next(rows)
        indices=[header.index(c) for c in columns]
        result=[[r[i].strip() if isinstance(r[i],str) else r[i] for i in indices]
                for r in rows if r[indices[0]] is not None]
    finally:book.close()
    sources.append(dict(file=relative,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return result
def save(name,columns,rows):
    with (OUT/name).open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f);writer.writerow(columns);writer.writerows(rows)
weather=read("MonthData/weatherresults.xlsx",["month","freezedays","hotdays","prcp","tavg","tmax","tmin","humidity"])
crime=dict(read("MonthData/CrimeByMonth2007.xlsx",["month","crimeindex"]))
months={name:i for i,name in enumerate(calendar.month_name) if i}
month_rows=sorted([[months.get(r[0],r[0])]+r[1:]+[crime[r[0]]] for r in weather],key=lambda r:int(r[0]))
assert len(month_rows)==12 and {int(r[0]) for r in month_rows}==set(range(1,13))
save("month_context.csv",["DMONTH","freezedays","hotdays","prcp","tavg","tmax","tmin","humidity","crimeindex"],month_rows)
crashes=dict(read("DayData/carcrashes.xlsx",["day","crashperc"]))
substance=dict(read("DayData/substancedeaths.xlsx",["day","SubstanceDeathsPer100k"]))
happiness=dict(read("DayData/happiness.xlsx",["day","happiness"]))
assert set(crashes)==set(substance)==set(happiness)==set(calendar.day_name)
save("weekday_context.csv",["DDAY","crashperc","substancedeaths","happiness"],
     [[d,crashes[d],substance[d],happiness[d]] for d in sorted(crashes)])
save("temporal_context_sources.csv",["file","sha256"],[[r["file"],r["sha256"]] for r in sources])
print("Prepared context for 12 months and seven weekdays.")
