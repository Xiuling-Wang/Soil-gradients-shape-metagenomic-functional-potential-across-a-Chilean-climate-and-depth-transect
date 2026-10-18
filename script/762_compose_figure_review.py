"""Assemble complete review figures from unchanged or newly rendered vector panels."""
from pathlib import Path
import os
import copy
import io
import json
import tempfile
import numpy as np
from pypdf import PdfReader, PdfWriter, PageObject, Transformation
from pypdf.generic import RectangleObject
import pypdfium2 as pdfium
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT = Path(os.environ.get("CHILE_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
AUDIT = ROOT / "analysis/figure_table_refinement_2026-10-18"
OUT = ROOT / "Environmental_Microbiome/Figure_table_review_2026-10-18"
NATIVE = AUDIT / "native"
BEFORE = AUDIT / "before"
MM = 72/25.4
PREVIEW = Path("/private/tmp/chile_figure_review_20261018")


def page(path, number=0):
    return PdfReader(path).pages[number]


def crop_top(p, top, bottom, left=0, right=None):
    q = copy.deepcopy(p)
    right = right if right is not None else float(q.mediabox.width)
    h = float(q.mediabox.height)
    q.cropbox = RectangleObject([left,h-bottom,right,h-top])
    target = PageObject.create_blank_page(width=right-left,height=bottom-top)
    target.merge_transformed_page(q,Transformation().translate(-left,-(h-bottom)))
    return target


def overlay(w,h,draw):
    b=io.BytesIO();c=canvas.Canvas(b,pagesize=(w,h));draw(c);c.save();b.seek(0)
    return PdfReader(b).pages[0]


def write(p,path):
    writer=PdfWriter();writer.add_page(p);writer.write(path)


def ink(p):
    b=io.BytesIO();w=PdfWriter();w.add_page(p);w.write(b)
    doc=pdfium.PdfDocument(b.getvalue());pg=doc[0]
    image=np.asarray(pg.render(scale=2).to_pil().convert("RGB"))
    mask=(image.min(axis=2)<248)
    ys,xs=np.where(mask)
    width,height=float(p.mediabox.width),float(p.mediabox.height)
    return [float(xs.min()/image.shape[1]*width),float((image.shape[0]-ys.max()-1)/image.shape[0]*height),float((xs.max()+1)/image.shape[1]*width),float((image.shape[0]-ys.min())/image.shape[0]*height)]


def compose(name, width, parts, gaps=(), labels=None, top_padding=0):
    height=max(y+float(p.mediabox.height)*s for _,p,x,y,s in parts)+top_padding
    result=PageObject.create_blank_page(width=width,height=height)
    records=[]
    for ident,p,x,y,s in parts:
        sw,sh=float(p.mediabox.width),float(p.mediabox.height)
        result.merge_transformed_page(p,Transformation().scale(s).translate(x,y))
        box=ink(p)
        records.append({"id":ident,"role":"legend" if "legend" in ident else "panel","source_width_pt":sw,"source_height_pt":sh,
            "placed_box_pt":[x,y,x+sw*s,y+sh*s],"ink_box_pt":[x+box[0]*s,y+box[1]*s,x+box[2]*s,y+box[3]*s]})
    if labels:
        result.merge_page(overlay(width,height,lambda c:labels(c,width,height)))
    write(result,OUT/name)
    (AUDIT/f"geometry_{name[:-4]}.json").write_text(json.dumps({
        "figure_id":name[:-4],"aspect_tolerance_pct":.05,"source_panels_locked_before_composition":True,
        "uses_outer_boxes_for_alignment":False,"panels":records,"gaps":list(gaps)},indent=2))


def main():
    PREVIEW.mkdir(exist_ok=True)
    pdfmetrics.registerFont(TTFont("Arial", str(Path(os.environ.get("CHILE_FONT_DIR", "/System/Library/Fonts/Supplemental"))/"Arial.ttf")))
    pdfmetrics.registerFont(TTFont("Arial-Bold", str(Path(os.environ.get("CHILE_FONT_DIR", "/System/Library/Fonts/Supplemental"))/"Arial Bold.ttf")))
    p=page(BEFORE/"Figure_1.pdf");w,h=float(p.mediabox.width),float(p.mediabox.height)
    def scale_label(c):
        c.setFont("Arial",8.2);c.setFillColorRGB(.2,.23,.25);c.drawCentredString(374,h-16,"Square-root scale")
    p.merge_page(overlay(w,h,scale_label));write(p,OUT/"Figure_1.pdf")

    p=page(BEFORE/"Figure_2.pdf")
    a=crop_top(p,15,275);leg=crop_top(p,281,313);c=crop_top(p,320,445)
    compose("Figure_2.pdf",518,[("C",c,0,0,1),("legend",leg,0,130,1),("AB",a,0,157,1)],
        [{"panel_a":"AB","panel_b":"C","axis":"y","kind":"vertical","min_pt":25,"max_pt":65}])

    a=crop_top(page(BEFORE/"Figure_5.pdf"),0,213)
    b=page(NATIVE/"Figure_5_B.pdf");width=float(b.mediabox.width)
    compose("Figure_5.pdf",width,[("B",b,0,0,1),("A",a,0,float(b.mediabox.height)+7,width/float(a.mediabox.width))],
        [{"panel_a":"A","panel_b":"B","axis":"y","kind":"vertical","min_pt":6,"max_pt":25}])

    a=page(AUDIT/"figure6/native_panels/Figure_6_A.pdf")
    b=page(AUDIT/"figure6/native_panels/Figure_6_B.pdf")
    width=float(a.mediabox.width)
    padded=[]
    for letter,p in [("A",a),("B",b)]:
        pw,ph=float(p.mediabox.width),float(p.mediabox.height)
        def erase_tag(c):
            c.setFillColorRGB(1,1,1);c.rect(0,ph-20,16,20,fill=1,stroke=0)
        p.merge_page(overlay(pw,ph,erase_tag))
        q=PageObject.create_blank_page(width=pw+12,height=ph+12)
        q.merge_transformed_page(p,Transformation().translate(6,6))
        def draw_tag(c,letter=letter):
            c.setFont("Arial-Bold",15.5);c.drawString(6,ph+12-17,letter)
        q.merge_page(overlay(pw+12,ph+12,draw_tag))
        write(q,NATIVE/f"Figure_6_{letter}_padded.pdf")
        padded.append(q)
    a,b=padded;s=width/float(a.mediabox.width)
    y_a=ink(b)[3]*s+26-ink(a)[1]*s
    compose("Figure_6.pdf",width,[("B",b,0,0,s),("A",a,0,y_a,s)],
        [{"panel_a":"A","panel_b":"B","axis":"y","kind":"vertical","min_pt":22,"max_pt":30}])

    p=page(BEFORE/"Figure_S1.pdf");w,h=float(p.mediabox.width),float(p.mediabox.height)
    bottom=2.65*72/2.54
    def erase_c(c):
        c.setFillColorRGB(1,1,1);c.rect(0,bottom,264,249,fill=1,stroke=0)
    p.merge_page(overlay(w,h,erase_c))
    p.merge_transformed_page(page(NATIVE/"Figure_S1_C.pdf"),Transformation().translate(0,bottom))
    write(p,OUT/"Figure_S1.pdf")

    for letter in ["A","B"]:
        plot=page(NATIVE/f"Figure_S2_{letter}.pdf")
        if letter=="A":
            assert ink(plot)[1]>37
            plot=crop_top(plot,0,float(plot.mediabox.height)-35)
        legend=page(NATIVE/f"Figure_S2_{letter}_legend.pdf")
        width=float(plot.mediabox.width)
        def labels(c,w,h,letter=letter):
            c.setFont("Arial-Bold",15.5);c.drawString(3,h-15,letter)
            if letter=="B":
                c.setFont("Arial",9);c.drawRightString(w-4,h-11,"Fig. S2 (continued)")
        compose(f"Figure_S2_{letter}.pdf",width,[("legend",legend,0,0,1),(letter,plot,0,float(legend.mediabox.height)+6,1)],
            gaps=[{"panel_a":letter,"panel_b":"legend","axis":"y","kind":"vertical","min_pt":12,"max_pt":36}],
            labels=labels,top_padding=18 if letter=="B" else 0)
    writer=PdfWriter()
    for letter in ["A","B"]:writer.add_page(page(OUT/f"Figure_S2_{letter}.pdf"))
    writer.write(OUT/"Figure_S2.pdf")
    # Separate panels remain internal; the delivered S2 is always the complete figure.
    for letter in ["A","B"]:
        (OUT/f"Figure_S2_{letter}.pdf").replace(NATIVE/f"Figure_S2_{letter}_composed.pdf")

    for path in sorted(OUT.glob("*.pdf")):
        doc=pdfium.PdfDocument(path)
        for i in range(len(doc)):
            img=doc[i].render(scale=2.1).to_pil()
            img.save(PREVIEW/f"{path.stem}_{i+1}.png")
    print(json.dumps({"output":str(OUT),"preview":str(PREVIEW),"figures":[p.name for p in sorted(OUT.glob('*.pdf'))]}))


if __name__=="__main__":main()
