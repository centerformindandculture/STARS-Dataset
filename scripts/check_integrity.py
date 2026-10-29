"""Create or verify SHA-256 checksums for public analysis files."""
from pathlib import Path
import argparse, hashlib
ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/"docs/CHECKSUMS.sha256"
def files():
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file() or p==MANIFEST:continue
        rel=p.relative_to(ROOT)
        if set(rel.parts)&{".git",".venv",".Rproj.user","__pycache__"}:continue
        if p.suffix.lower() in {".docx",".pdf",".pyc"}:continue
        if p.suffix.lower()==".xlsx" and rel.parts[:2]!=("data","external"):continue
        if rel.parts[:2] in {("data","raw"),("outputs","local"),("outputs","reproduced")}:continue
        if p.name in {".Rhistory",".RData",".DS_Store","Thumbs.db","desktop.ini"}:continue
        yield p
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true",help="Generate release checksums")
    args=parser.parse_args()
    actual={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files()}
    if args.write:
        MANIFEST.write_text("".join(f"{v}  {k}\n" for k,v in actual.items()),encoding="utf-8")
        print(f"Wrote checksums for {len(actual)} published files.");return
    expected={}
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        value,key=line.split("  ",1);expected[key]=value
    failures=[f"Missing or changed: {k}" for k,v in expected.items() if actual.get(k)!=v]
    failures += [f"Not in manifest: {k}" for k in actual.keys()-expected.keys()]
    if failures:raise SystemExit("\n".join(failures))
    print(f"Verified {len(expected)} published files.")
if __name__=="__main__":main()
