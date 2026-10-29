from pathlib import Path
import os, json, math, copy, re, xml.etree.ElementTree as ET
from PIL import Image, ImageDraw, ImageFont
root=Path(__file__).resolve().parents[1]
out=root/Path(os.environ.get('STARS_OUTPUT_DIR','outputs'))/'clustering/geographic'
NS='http://www.w3.org/2000/svg';ET.register_namespace('',NS)
original=ET.parse(out/'county_k3_map.svg').getroot()
paths=original.findall('.//{'+NS+'}path')
assert len(paths)==56
W,H=1800,1000
svg=ET.Element('{'+NS+'}svg',{'viewBox':f'0 0 {W} {H}','role':'img','aria-label':'Montana counties colored by three contextual clusters'})
ET.SubElement(svg,'{'+NS+'}rect',{'width':'100%','height':'100%','fill':'white'})
map_group=ET.SubElement(svg,'{'+NS+'}g',{'transform':'translate(0,-110)'})
for path in paths:map_group.append(copy.deepcopy(path))
colors={1:'#E69F00',2:'#56B4E9',3:'#009E73'};counts={1:28,2:6,3:22}
legend=ET.SubElement(svg,'{'+NS+'}g',{'font-family':'Arial, sans-serif','font-size':'38','fill':'#172B36'})
for g in (1,2,3):
 x=100+(g-1)*560
 ET.SubElement(legend,'{'+NS+'}rect',{'x':str(x),'y':'935','width':'40','height':'36','fill':colors[g]})
 text=ET.SubElement(legend,'{'+NS+'}text',{'x':str(x+56),'y':'966'})
 text.text=f'Cluster {g} ({counts[g]} counties)'
ET.ElementTree(svg).write(out/'county_k3_map_legend_only.svg',encoding='utf-8',xml_declaration=True)
assert [p.attrib['d'] for p in paths]==[p.attrib['d'] for p in map_group]
assert len(svg.findall('.//{'+NS+'}text'))==3
# Render the identical map paths at twice native resolution for publication.
S=2;im=Image.new('RGB',(W*S,H*S),'white');dr=ImageDraw.Draw(im)
for path in paths:
 rings=[]
 for subpath in path.attrib['d'].split('M ')[1:]:
  pts=[tuple(map(float,p)) for p in re.findall(r'([-\d.]+),([-\d.]+)',subpath)]
  rings.append([(x*S,(y-110)*S) for x,y in pts])
 mask=Image.new('L',im.size);md=ImageDraw.Draw(mask)
 # The source uses even-odd filling; toggle each ring, retaining island holes.
 from PIL import ImageChops
 for ring in rings:
  one=Image.new('1',im.size);ImageDraw.Draw(one).polygon(ring,fill=1)
  mask=ImageChops.difference(mask,one.convert('L'))
 im.paste(path.attrib['fill'],mask=mask)
 for ring in rings:dr.line(ring+[ring[0]],fill='white',width=round(1.8*S),joint='curve')
font=None
for name in ['C:/Windows/Fonts/arial.ttf','arial.ttf','DejaVuSans.ttf']:
 try:font=ImageFont.truetype(name,38*S);break
 except OSError:pass
assert font is not None
for g in (1,2,3):
 x=100+(g-1)*560;label=f'Cluster {g} ({counts[g]} counties)'
 assert x+56+dr.textlength(label,font=font)/S<W-40
 dr.rectangle((x*S,935*S,(x+40)*S,971*S),fill=colors[g])
 dr.text(((x+56)*S,966*S),label,fill='#172B36',font=font,anchor='ls')
im.save(out/'county_k3_map_legend_only.png',dpi=(300,300))
print('Created legend-only SVG and 3600 x 2000 PNG with 56 counties and 38px legend labels.')
