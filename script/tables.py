"""Apply approved table presentation edits without changing observations."""
from pathlib import Path
from zipfile import ZipFile
from copy import deepcopy
import argparse
from lxml import etree as E

URI='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
NS={'s':URI}
Q=lambda name:'{'+URI+'}'+name
WIDTHS=dict(zip('GHIJKLMNOPQRSTUVW',[16,12,15,10,17,13,15,15,15,15,24,24,24,12,13,12,12]))
FORMATS={'H':'0.0','J':'0.00','K':'0.00','L':'0.00','Q':'0.00','R':'0.00','S':'0.00'}
HEADERS={'B3':'Completeness (%)','C3':'Contamination (%)','D3':'Genome size (bp)',
         'E3':'GC content (fraction)','F3':'N50 (bp)'}

def serialize(root):
    return E.tostring(root,encoding='UTF-8',xml_declaration=True,standalone=True)

def revise(source,dest,kind):
    if dest.exists(): raise ValueError('Output already exists: '+str(dest))
    with ZipFile(source) as z:
        infos=z.infolist();original={n:z.read(n) for n in z.namelist()}
    result=dict(original)
    if kind=='tables':
        styles=E.fromstring(original['xl/styles.xml'])
        xfs=styles.find('s:cellXfs',NS);formats=styles.find('s:numFmts',NS)
        known={n.get('formatCode'):int(n.get('numFmtId')) for n in formats}
        for fmt in dict.fromkeys(FORMATS.values()):
            if fmt not in known:
                idx=max([163]+list(known.values()))+1
                E.SubElement(formats,Q('numFmt'),numFmtId=str(idx),formatCode=fmt);known[fmt]=idx
        formats.set('count',str(len(formats)))
        s1=E.fromstring(original['xl/worksheets/sheet2.xml'])
        cells={c.get('r'):c for c in s1.findall('.//s:sheetData/s:row/s:c',NS)}
        for letter,width in WIDTHS.items():
            idx=ord(letter)-64
            col=next(c for c in s1.findall('s:cols/s:col',NS) if int(c.get('min'))<=idx<=int(c.get('max')))
            assert col.get('min')==col.get('max')
            col.set('width',str(width));col.set('customWidth','1')
        cache={}
        for letter,fmt in FORMATS.items():
            for row in range(4,21):
                c=cells[f'{letter}{row}'];old=int(c.get('s','0'));key=(old,known[fmt])
                if key not in cache:
                    xf=deepcopy(xfs[old]);xf.set('numFmtId',str(known[fmt]));xf.set('applyNumberFormat','1')
                    cache[key]=len(xfs);xfs.append(xf)
                c.set('s',str(cache[key]))
        xfs.set('count',str(len(xfs)));result['xl/styles.xml']=serialize(styles)
        s4=E.fromstring(original['xl/worksheets/sheet8.xml'])
        for sheet,height,path in [(s1,66,'xl/worksheets/sheet2.xml'),(s4,56,'xl/worksheets/sheet8.xml')]:
            row=sheet.find("s:sheetData/s:row[@r='3']",NS);row.set('ht',str(height));row.set('customHeight','1')
            if sheet is s4:
                cs={c.get('r'):c for c in sheet.findall('.//s:sheetData/s:row/s:c',NS)}
                for addr,text in HEADERS.items():
                    c=cs[addr]
                    for child in list(c):
                        if child.tag in (Q('v'),Q('is')): c.remove(child)
                    c.set('t','inlineStr');E.SubElement(E.SubElement(c,Q('is')),Q('t')).text=text
            result[path]=serialize(sheet)
    for path,data in original.items():
        if not(path.startswith('xl/tables/') and path.endswith('.xml')): continue
        tab=E.fromstring(data)
        if kind=='data':
            tab.find('s:tableStyleInfo',NS).set('name','TableStyleLight9');result[path]=serialize(tab)
        elif tab.get('name')=='V40_TableS4':
            cols=tab.find('s:tableColumns',NS)
            for addr,header in HEADERS.items(): cols[ord(addr[0])-65].set('name',header)
            result[path]=serialize(tab)
    dest.parent.mkdir(parents=True,exist_ok=True)
    with ZipFile(dest,'w') as z:
        for info in infos: z.writestr(info,result[info.filename])

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,required=True,help='Directory containing the two original workbooks')
    p.add_argument('--out',type=Path,required=True,help='New output directory')
    a=p.parse_args()
    if a.out.exists(): p.error('--out must be a new directory')
    names=[('Tables_S1-S4.xlsx','tables'),('Supplementary_Data_S1.xlsx','data')]
    for name,_ in names:
        if not (a.inputs/name).is_file(): p.error('Missing '+name)
    for name,kind in names: revise(a.inputs/name,a.out/name,kind)

if __name__=='__main__': main()
