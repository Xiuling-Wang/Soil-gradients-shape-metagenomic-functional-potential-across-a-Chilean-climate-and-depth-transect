"""Graphical synthesis of the same frozen results and selected MAG combinations."""
from pathlib import Path
import os
import csv
import json
from reportlab.graphics import renderPDF, renderSVG
from reportlab.graphics.shapes import Drawing, Line, Rect, String, Circle, Wedge
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor

ROOT = Path(os.environ.get("CHILE_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
OUT = ROOT / "Environmental_Microbiome/Figure_table_review_2026-10-18_round2"
AUDIT = ROOT / "analysis/figure_table_refinement_2026-10-18_round2"
PARENT = ROOT / "Environmental_Microbiome/Coauthor_Review_2026-10-01_XW"
W, H = 510, 461
INK, MUTED, RULE = "#22282D", "#526068", "#B7C3C8"
BLUE, TEAL, PINK, GREEN, GOLD = "#3276A8", "#238491", "#B44583", "#468749", "#947328"


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    for name, filename in [("Arial", "Arial.ttf"), ("Arial-Bold", "Arial Bold.ttf"), ("Arial-Italic", "Arial Italic.ttf")]:
        pdfmetrics.registerFont(TTFont(name, str(Path(os.environ.get("CHILE_FONT_DIR", "/System/Library/Fonts/Supplemental"))/filename)))
    d, boxes, lines = Drawing(W, H), [], []

    def text(x, y, value, size=9.3, font="Arial", colour=INK, anchor="start"):
        width = pdfmetrics.stringWidth(value, font, size)
        left = x - (width / 2 if anchor == "middle" else width if anchor == "end" else 0)
        box = [left, y-size*.8, left+width, y+size*.22]
        assert left >= 0 and left+width <= W-3 and box[1] >= 0 and box[3] <= H-2, (value, box)
        boxes.append({"text": value, "box_pt": box, "font_pt": size, "font": font})
        d.add(String(x, H-y, value, fontName=font, fontSize=size, fillColor=HexColor(colour), textAnchor=anchor))
        return width

    def rich(x, y, parts, size=9.1):
        for value, font in parts:
            x += text(x, y, value, size, font)

    def line(x0, y0, x1, y1, colour=RULE, width=.8):
        d.add(Line(x0, H-y0, x1, H-y1, strokeColor=HexColor(colour), strokeWidth=width))
        lines.append([x0,y0,x1,y1])

    def rect(x, y, width, height, fill, stroke=None, radius=0):
        d.add(Rect(x,H-y-height,width,height,rx=radius,ry=radius,fillColor=HexColor(fill),
                   strokeColor=HexColor(stroke) if stroke else None,strokeWidth=.7))

    def gene(x, y, width, label, colour=GREEN, fill="#F1F7F1", size=9.7):
        rect(x,y,width,21,fill,colour,3)
        text(x+width/2,y+14,label,size,"Arial-Italic",anchor="middle")

    def module(x, y, genes, function, colour, second=None):
        line(x,y,x+7,y,colour)
        d.add(Circle(x,H-y,2.4,fillColor=HexColor(colour),strokeColor=None))
        text(x+12,y+3,genes,9.5,"Arial-Italic" if genes[0].islower() else "Arial")
        if second:
            text(x+12,y+17,second,9.5,"Arial-Italic")
        if function:
            text(x+12,y+(31 if second else 17),function,8.3,colour=MUTED)

    rows = list(csv.DictReader((PARENT/"Supplementary_Data/Category_pit_endpoints.tsv").open(),delimiter="\t"))
    groups = {}
    for row in rows:
        groups.setdefault(row["category"],[]).append(row)
    declining = [key for key,values in groups.items() if len(values)==3 and all(float(r["deep_TPM"])<float(r["surface_TPM"]) for r in values)]
    assert len(groups)==16 and len(declining)==6

    text(2,16,"A",14,"Arial-Bold")
    text(24,16,"Regional differentiation and selective depth patterns",10.6,"Arial-Bold")
    for x,title,r2,p in [(24,"KO composition: four surface sites",.479,"0.0375"),(281,"KO composition: NA depth profiles",.389,"0.361")]:
        text(x,40,title,10.4,"Arial-Bold")
        # The same existing R-squared is encoded on a fixed 0-1 length scale.
        rect(x,50,125,5,"#E8EDF0")
        rect(x,50,125*r2,5,BLUE)
        rich(x,76,[("R","Arial-Italic"),(f"\u00b2 = {r2:.3f}   Adjusted ","Arial"),("p","Arial-Italic"),(f" = {p}","Arial")])
    text(24,98,"Higher surface-site mean category TPM",8.8,colour=MUTED)
    for y,value,colour in [(114,"Hydrogenases: AZ, SG",TEAL),(130,"Putative CO oxidation: LC, NA",MUTED),(146,"Sulfur oxidation: NA",GOLD)]:
        d.add(Circle(27,H-y+3,2.5,fillColor=HexColor(colour),strokeColor=None))
        text(35,y,value,9.0)

    d.add(Circle(307,H-119,23,fillColor=HexColor("#E8EDF0"),strokeColor=None))
    d.add(Wedge(307,H-119,23,90,90+360*6/16,fillColor=HexColor(TEAL),strokeColor=None))
    d.add(Circle(307,H-119,17,fillColor=HexColor("#FFFFFF"),strokeColor=None))
    text(307,122,"6/16",10,"Arial-Bold",anchor="middle")
    text(342,107,"Categories lower in",9.2,"Arial-Bold")
    text(342,121,"every soil profile",9.2,"Arial-Bold")
    text(342,139,"5\u201310 to 60\u201380 cm",9,colour=MUTED)
    text(281,156,"N, energy and weathering categories",8.3,colour=MUTED)

    text(2,184,"B",14,"Arial-Bold")
    text(24,184,"Shared CBB genes, distinct genomic combinations",10.6,"Arial-Bold")
    text(255,205,"CBB-associated neighbourhood in all three MAGs",9.2,"Arial-Bold",anchor="middle")
    line(179,226,331,226,GREEN,1)
    gene(183,215,43,"prk")
    gene(234,215,43,"rbcL")
    gene(285,215,43,"rbcS")
    # Undirected inclusion links, not reaction arrows or inferred gene orientation.
    line(255,237,255,248)
    line(84,248,426,248)
    for x in [84,255,426]:
        line(x,248,x,263)

    for x,mag,taxon,fill,colour,italic in [
        (84,"MAG122","Bradyrhizobium","#F4F9F0","#789857",True),
        (255,"MAG138","Actinophytocola","#F0F7FA","#6B99AC",True),
        (426,"MAG156","CF-13 (Dormibacterota)","#EFF5FA","#3276A8",False),
    ]:
        rect(x-72,263,144,40,fill,colour,5)
        text(x,279,mag,10.8,"Arial-Bold",anchor="middle")
        text(x,294,taxon,8.6,"Arial-Italic" if italic else "Arial",anchor="middle")

    line(84,303,84,312);line(18,312,84,312);line(18,312,18,413)
    module(18,326,"Sox genes","Sulfur oxidation",GOLD)
    module(18,366,"ureA / ureB / ureC","Urea hydrolysis",BLUE)
    module(18,413,"pqqB / pqqC","Organic-acid weathering",PINK,second="pqqD / pqqE")

    line(255,303,255,312);line(189,312,255,312);line(189,312,189,413)
    module(189,326,"Hydrogenase","homologues",TEAL)
    module(189,366,"nirB / nirD",None,BLUE,second="ureA / ureB / ureC")
    module(189,413,"pqqB / pqqC","Organic-acid weathering",PINK,second="pqqD / pqqE")

    line(426,303,426,330)
    text(426,344,"Local neighbourhood",9,colour=MUTED,anchor="middle")
    line(355,365,496,365,PINK,1)
    gene(360,354,55,"ppk1",PINK,"#FCF2F7")
    gene(438,354,55,"pstS",PINK,"#FCF2F7")
    line(387.5,376,387.5,391,PINK,.75)
    line(465.5,376,465.5,391,PINK,.75)
    for x,words in [(387.5,["Polyphosphate","synthesis"]),(465.5,["Phosphate","uptake"])]:
        for i,value in enumerate(words):
            text(x,405+i*13,value,8.5,anchor="middle")

    overlaps=[]
    for i,a in enumerate(boxes):
        ax,ay,ar,ab=a["box_pt"]
        for b in boxes[i+1:]:
            bx,by,br,bb=b["box_pt"]
            if min(ar,br)-max(ax,bx)>0.02 and min(ab,bb)-max(ay,by)>0.02:
                overlaps.append((a["text"],b["text"]))
    assert not overlaps,overlaps
    renderPDF.drawToFile(d,str(OUT/"Figure_7.pdf"))
    renderSVG.drawToFile(d,str(AUDIT/"Figure_7.editable.svg"))
    (AUDIT/"Figure7_QA.json").write_text(json.dumps({
        "native_size_pt":[W,H],"text_boxes":boxes,"lines":lines,"overlaps":overlaps,
        "italic_p_count":sum(b["text"]=="p" and b["font"]=="Arial-Italic" for b in boxes),
        "R2_bar_scale":[0,1],"statistical_results_unchanged":[.479,.389,.0375,.361],
        "categories":16,"recurring_declines":declining,"gene_orientation_not_inferred":True,
        "undirected_connections":"Shared neighbourhood membership and selected co-occurring traits, not metabolic flux",
        "new_analysis":False,"author_acceptance":"approved selection",
    },indent=2))
    print("Figure 7: vector diagram, same results, italic p, no text overlaps")


if __name__=="__main__":
    main()
