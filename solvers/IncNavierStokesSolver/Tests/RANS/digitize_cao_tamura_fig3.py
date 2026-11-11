#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Extract hollow rough-case vector markers; no PDF or invented profile in Git."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

PDF_SHA='adffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df'
# Labeled ticks in original PDF points (bottom-origin). Never infer a zero
# from an unlabeled plot border: shear panel's lower border is above y/delta=0.
PANELS={
 'a':dict(box=[100,240,450,615],legend=[129.6305,575.863],x_ticks=[[106.554,0],[222.474,1]],
          y_ticks=[[457.518,0],[500.562,100]],value='U/U_infinity',height='y_mm'),
 'b':dict(box=[265,405,450,615],legend=[320.0905,577.5245],x_ticks=[[273.171,0],[299.466,.001]],
          y_ticks=[[519.912,.5],[584.506,1]],value='-R_xz/U_infinity^2',height='y/delta'),
 'c':dict(box=[100,240,275,430],legend=[157.117,399.451],x_ticks=[[106.019,0],[233.999,3]],
          y_ticks=[[331.884,.5],[390.017,1]],value='sigma_u/u_tau',height='(y-d)/delta'),
 'd':dict(box=[265,405,275,430],legend=[323.1015,399.746],x_ticks=[[272.324,0],[335.926,1]],
          y_ticks=[[332.553,.5],[390.365,1]],value='sigma_vertical/u_tau',height='(y-d)/delta')}

def linear(x,ticks):
    (a,u),(b,v)=ticks;return u+(x-a)*(v-u)/(b-a)

def extract(pdf):
    import logging
    logging.getLogger('pypdf').setLevel(logging.ERROR)
    from pypdf import PdfReader
    from pypdf.generic import ContentStream
    assert hashlib.sha256(pdf.read_bytes()).hexdigest()==PDF_SHA
    reader=PdfReader(pdf);page=reader.pages[6];fill=None;path=[];markers=[]
    for values,op in ContentStream(page.get_contents(),reader).operations:
        if op==b'scn':fill=[float(x) for x in values]
        if op in (b'm',b'l',b'c',b'h'):path.append((op,list(map(float,values))))
        if op in (b'f*',b'f',b'S',b's',b'B',b'B*',b'n'):
            if len(path)==5 and path[0][0]==b'm' and all(p==b'c' for p,_ in path[1:]) and fill==[1.,1.,1.] and op in (b'f*',b'f'):
                pts=[path[0][1]]+[v[-2:] for _,v in path[1:]]
                lo=[min(p[i] for p in pts) for i in range(2)];hi=[max(p[i] for p in pts) for i in range(2)]
                if .8<hi[0]-lo[0]<5 and .8<hi[1]-lo[1]<5:markers.append([(a+b)/2 for a,b in zip(lo,hi)])
            path=[]
    rows=[]
    for panel,cal in PANELS.items():
        a,b,c,d=cal['box'];points=[p for p in markers if a<p[0]<b and c<p[1]<d
            and sum((p[i]-cal['legend'][i])**2 for i in range(2))>.25]
        assert len(points)==15,(panel,len(points))
        for i,(x,y) in enumerate(sorted(points,key=lambda p:p[1])):
            v=linear(x,cal['x_ticks']);height=linear(y,cal['y_ticks'])
            physical=height if panel=='a' else height*250+(3.5 if panel in ('c','d') else 0)
            rows.append(dict(panel=panel,point=i+1,pdf_x=x,pdf_y=y,value=v,value_normalization=cal['value'],
                raw_height=height,height_origin=cal['height'],y_mm=physical,y_over_H=physical/40,
                r_over_H=(physical-3.5)/40,kind='digitized hollow rough-case circle',page=7,figure='3('+panel+')'))
    return rows,float(page.mediabox.width),float(page.mediabox.height)

def main():
    p=argparse.ArgumentParser();p.add_argument('pdf',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--page-png',type=Path,help='optional independent page7 raster for source-panel overlays')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);rows,w,h=extract(a.pdf)
    with (a.output/'rough-fig3.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    resolution={}
    for name,cal in PANELS.items():
        dx=abs(linear(.15,cal['x_ticks'])-linear(0,cal['x_ticks']))
        dy=abs(linear(.15,cal['y_ticks'])-linear(0,cal['y_ticks']))*(1 if name=='a' else 250)
        resolution[name]=dict(value_reading_resolution=dx,height_mm_reading_resolution=dy)
    (a.output/'digitization.json').write_text(json.dumps(dict(schema='CaoTamuraFig3RoughV1',pdf_sha256=PDF_SHA,
        pdf_page=7,page_size_points=[w,h],panels=PANELS,resolution=resolution,
        source_selection='hollow white-filled four-Bezier circles, excluding legends; excludes smooth filled circles and Raupach crosses',
        reading_assumption='conservative +/-0.15 PDF-point center/tick reading resolution; not an experimental uncertainty or confidence interval',
        physical_conversion='delta=250mm and d=3.5mm assumed as reported; preserve panel-specific raw heights',
        data_plane='virtual crest with hill absent,5m downstream of tunnel contraction; not arbitrary upstream inlet',
        excluded_error_claim='probe voltage calibration +/-1.5% is not a universal velocity uncertainty'),indent=2)+'\n')
    if a.page_png:
        from PIL import Image,ImageDraw
        image=Image.open(a.page_png).convert('RGB');sx=image.width/w;sy=image.height/h
        panels=[]
        for name,cal in PANELS.items():
            lo,hi,bot,top=cal['box'];left=int((lo-15)*sx);upper=int((h-top-5)*sy)
            crop=image.crop((left,upper,int((hi+8)*sx),int((h-bot+16)*sy)))
            draw=ImageDraw.Draw(crop)
            for row in rows:
                if row['panel']!=name:continue
                x=row['pdf_x']*sx-left;y=(h-row['pdf_y'])*sy-upper
                draw.line((x-3,y,x+3,y),fill='red',width=1);draw.line((x,y-3,x,y+3),fill='red',width=1)
            panels.append(crop)
        width=max(p.width for p in panels);height=max(p.height for p in panels)
        canvas=Image.new('RGB',(2*width,2*height),'white')
        for i,panel in enumerate(panels):canvas.paste(panel,((i%2)*width,(i//2)*height))
        canvas.save(a.output/'source-overlay.png')

if __name__=='__main__':main()
