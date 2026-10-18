"""Render the approved figure selection in a new, isolated workspace."""
from pathlib import Path
import argparse
import importlib.util
import os
import shutil
import subprocess
import sys

PARENT = 'Environmental_Microbiome/Coauthor_Review_2026-10-01_XW'
A1 = 'analysis/figure_table_refinement_2026-10-18'
A2 = 'analysis/figure_table_refinement_2026-10-18_round2'
D1 = 'Environmental_Microbiome/Figure_table_review_2026-10-18'
D2 = 'Environmental_Microbiome/Figure_table_review_2026-10-18_round2'
INPUTS = [
    'data/metadata/env_2024_meta.tsv',
    'data/processed/metagenome_output/MAGs',
    PARENT+'/Main_Figures_PDF', PARENT+'/Supplementary_Figures_PDF',
    PARENT+'/Supplementary_Data',
    'Environmental_Microbiome/V35/Supplementary_Data/AmoA_sequences',
    'analysis/V42_ecological_story_revision_2026-08-26/source_data',
    'analysis/V53_author_comments_revision_2026-08-28/source_data',
    'analysis/V59_KO_ordination_architecture_2026-08-29/source_data',
    'analysis/V05_comment_review_2026-08-17/figure2_reanalysis',
    'analysis/exploratory_KO_differential_abundance_2026-08-29/KO_DA_README.md',
]

def load(root, name):
    spec = importlib.util.spec_from_file_location('figure_'+name[:3], root/'script'/name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def prepare(inputs, work, code):
    if work.exists():
        raise SystemExit('Choose a new --work directory; existing files are never overwritten.')
    missing = [p for p in INPUTS if not (inputs/p).exists()]
    if missing:
        raise SystemExit('Missing private/local inputs:\n'+'\n'.join(missing))
    work.mkdir(parents=True)
    shutil.copytree(code/'script', work/'script', ignore=shutil.ignore_patterns('__pycache__'))
    for rel in INPUTS:
        dest = work/rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.symlink_to(inputs/rel, target_is_directory=(inputs/rel).is_dir())
    for rel in [A1+'/native', A1+'/before', A2+'/native', D1, D2]:
        (work/rel).mkdir(parents=True, exist_ok=True)
    for group in ['Main_Figures_PDF', 'Supplementary_Figures_PDF']:
        for p in (inputs/PARENT/group).glob('*.pdf'):
            shutil.copy2(p, work/A1/'before'/p.name)

def compose(work):
    from pypdf import PdfWriter, PageObject, Transformation
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    m = load(work, '762_compose_figure_review.py')
    fontdir = Path(os.environ.get('CHILE_FONT_DIR', '/System/Library/Fonts/Supplemental'))
    for name, filename in [('Arial','Arial.ttf'), ('Arial-Bold','Arial Bold.ttf')]:
        pdfmetrics.registerFont(TTFont(name, str(fontdir/filename)))
    m.OUT = work/D2
    before, native = work/A1/'before', work/A1/'native'
    p = m.page(before/'Figure_1.pdf'); w,h = float(p.mediabox.width),float(p.mediabox.height)
    def label(c):
        c.setFont('Arial',8.2); c.setFillColorRGB(.2,.23,.25)
        c.drawCentredString(374,h-16,'Square-root scale')
    p.merge_page(m.overlay(w,h,label));m.write(p,m.OUT/'Figure_1.pdf')
    p = m.page(before/'Figure_2.pdf')
    a,leg,c = m.crop_top(p,15,275),m.crop_top(p,281,313),m.crop_top(p,320,445)
    m.compose('Figure_2.pdf',518,[('C',c,0,0,1),('legend',leg,0,130,1),('AB',a,0,157,1)])
    for n in [3,4]:
        shutil.copy2(work/D1/f'Figure_{n}.pdf',m.OUT/f'Figure_{n}.pdf')
    a = m.crop_top(m.page(before/'Figure_5.pdf'),0,210)
    b = m.crop_top(m.page(work/A2/'native/Figure_5_B_revised.pdf'),20,444)
    p = PageObject.create_blank_page(width=487,height=639)
    p.merge_transformed_page(a,Transformation().translate(0,429));p.merge_page(b)
    m.write(p,m.OUT/'Figure_5.pdf')
    panels=[]
    for letter in ['A','B']:
        p=m.page(work/A1/f'figure6/native_panels/Figure_6_{letter}.pdf')
        pw,ph=float(p.mediabox.width),float(p.mediabox.height)
        def erase(c):
            c.setFillColorRGB(1,1,1);c.rect(0,ph-20,16,20,fill=1,stroke=0)
        p.merge_page(m.overlay(pw,ph,erase))
        q=PageObject.create_blank_page(width=pw+12,height=ph+12)
        q.merge_transformed_page(p,Transformation().translate(6,6))
        def tag(c,letter=letter):
            c.setFont('Arial-Bold',15.5);c.drawString(6,ph+12-17,letter)
        q.merge_page(m.overlay(pw+12,ph+12,tag));panels.append(q)
    a,b=panels; width=float(a.mediabox.width)-12; s=width/float(a.mediabox.width)
    y=m.ink(b)[3]*s+26-m.ink(a)[1]*s
    m.compose('Figure_6.pdf',width,[('B',b,0,0,s),('A',a,0,y,s)])
    p=m.page(before/'Figure_S1.pdf');w,h=float(p.mediabox.width),float(p.mediabox.height)
    bottom=2.65*72/2.54
    def erase(c):
        c.setFillColorRGB(1,1,1);c.rect(0,bottom,264,249,fill=1,stroke=0)
    p.merge_page(m.overlay(w,h,erase))
    p.merge_transformed_page(m.page(native/'Figure_S1_C.pdf'),Transformation().translate(0,bottom))
    m.write(p,m.OUT/'Figure_S1.pdf')
    writer=PdfWriter()
    for letter in ['A','B']:
        writer.add_page(m.page(work/A2/f'native/Figure_S2_{letter}_page_revised.pdf'))
    writer.write(m.OUT/'Figure_S2.pdf')
    shutil.copy2(before/'Figure_S3.pdf',m.OUT/'Figure_S3.pdf')
    print('Complete figure selection:',m.OUT)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,required=True,help='Local input tree; no data are included in this repository')
    p.add_argument('--work',type=Path,required=True,help='New isolated output workspace')
    p.add_argument('--check',action='store_true',help='Check input paths only; write nothing')
    args=p.parse_args();inputs=args.inputs.resolve();work=args.work.resolve()
    if args.check:
        missing=[rel for rel in INPUTS if not (inputs/rel).exists()]
        for rel in missing: print('Missing:',rel)
        raise SystemExit(1 if missing else 0)
    code=Path(__file__).resolve().parents[1]
    prepare(inputs,work,code)
    os.environ['CHILE_PROJECT_ROOT']=str(work)
    def run(*args): subprocess.run(args,cwd=work,check=True)
    run('Rscript',str(work/'script/760_refine_native_figures.R'),'heatmaps','figureS1','figure6')
    run('Rscript',str(work/'script/766_author_figure_feedback_native.R'),'figure5')
    run('Rscript',str(work/'script/766_author_figure_feedback_native.R'),'figureS2')
    run(sys.executable,str(work/'script/767_figure7_gene_combination_diagram.py'))
    compose(work)

if __name__=='__main__': main()
