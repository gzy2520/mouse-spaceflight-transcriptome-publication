#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Apply precise tracked changes (WordprocessingML w:ins and w:del) to:
1. /Users/gzy2520/Desktop/mouse_Methods.docx
   - Italicize all mouse gene symbols in repair pathway modules (P18, 19, 20, 26, 27, 28, 29).
   - Clean up punctuation glitch in P36 (move '.' after parenthesis).
2. /Users/gzy2520/Downloads/Mouse_flight_with_legends_SB.docx
   - P12: Mre11A -> Mre11a
   - P33: was -> were (after Nhej1 and Paxx)
   - P33: increase -> increases (that increase -> that increases)
   - P33: ends -> end (ends ligation -> end ligation)
   - P38: occured -> occurred
   - P58: Ercc6 (plain) -> Ercc6 (italic) in Fig 5 legend
"""

import os
import re
import copy
import shutil
import zipfile
import tempfile
import lxml.etree as etree

NS = {
    'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
    'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math'
}

def repack_docx(source_dir, output_file):
    tmp_out = output_file + ".tmp.docx"
    with zipfile.ZipFile(tmp_out, 'w', zipfile.ZIP_DEFLATED) as z_out:
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                full_path = os.path.join(root, file)
                rel_path = os.path.relpath(full_path, source_dir)
                z_out.write(full_path, rel_path)
    shutil.move(tmp_out, output_file)

def make_del(text, rPr_base, cur_id, author="gzy", date="2026-09-10T12:45:00Z"):
    del_el = etree.Element(f"{{{NS['w']}}}del", {
        f"{{{NS['w']}}}id": str(cur_id),
        f"{{{NS['w']}}}author": author,
        f"{{{NS['w']}}}date": date
    })
    r = etree.SubElement(del_el, f"{{{NS['w']}}}r")
    if rPr_base is not None:
        r.append(copy.deepcopy(rPr_base))
    t = etree.SubElement(r, f"{{{NS['w']}}}delText")
    if text.startswith(' ') or text.endswith(' '):
        t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    t.text = text
    return del_el

def make_ins(text, rPr_base, cur_id, italic=False, author="gzy", date="2026-09-10T12:45:00Z"):
    ins_el = etree.Element(f"{{{NS['w']}}}ins", {
        f"{{{NS['w']}}}id": str(cur_id),
        f"{{{NS['w']}}}author": author,
        f"{{{NS['w']}}}date": date
    })
    r = etree.SubElement(ins_el, f"{{{NS['w']}}}r")
    rPr = copy.deepcopy(rPr_base) if rPr_base is not None else etree.Element(f"{{{NS['w']}}}rPr")
    if italic:
        if rPr.find('w:i', NS) is None:
            etree.SubElement(rPr, f"{{{NS['w']}}}i")
        if rPr.find('w:iCs', NS) is None:
            etree.SubElement(rPr, f"{{{NS['w']}}}iCs")
    r.append(rPr)
    t = etree.SubElement(r, f"{{{NS['w']}}}t")
    if text.startswith(' ') or text.endswith(' '):
        t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    t.text = text
    return ins_el

def make_plain_r(text, rPr_base):
    r = etree.Element(f"{{{NS['w']}}}r")
    if rPr_base is not None:
        r.append(copy.deepcopy(rPr_base))
    t = etree.SubElement(r, f"{{{NS['w']}}}t")
    if text.startswith(' ') or text.endswith(' '):
        t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    t.text = text
    return r

# ==============================================================================
# Process Document 1: mouse_Methods.docx
# ==============================================================================
def process_methods_doc(docx_path):
    print(f"Processing Methods doc: {docx_path}")
    temp_dir = tempfile.mkdtemp(prefix="methods_docx_")
    with zipfile.ZipFile(docx_path, 'r') as z:
        z.extractall(temp_dir)
    
    xml_path = os.path.join(temp_dir, 'word', 'document.xml')
    tree = etree.parse(xml_path)
    root = tree.getroot()
    body = root.find('w:body', NS)
    paras = body.findall('w:p', NS)
    
    # Determine max ID
    ids = []
    for el in root.iter():
        for k, v in el.attrib.items():
            if k.endswith('id'):
                try:
                    ids.append(int(v))
                except ValueError:
                    pass
    cur_id = max(ids) + 1 if ids else 500
    
    # 1. Italicize genes in P18, P19, P20, P26, P27, P28, P29
    genes = [
        'Nhej1', 'Paxx', 'Xrcc6', 'Xrcc5',
        'Brca1', 'Bard1', 'Blm', 'Rad51',
        'Parp1', 'Polq', 'Lig1', 'Lig3',
        'Ung', 'Ogg1', 'Neil1',
        'Xpc', 'Rad23b', 'Cetn2', 'Ddb1', 'Ddb2', 'Ercc6', 'Ercc8',
        'Msh2', 'Msh3',
        'Fancd2', 'Fanci'
    ]
    pattern = re.compile(r'\b(' + '|'.join(genes) + r')\b')
    
    target_indices = [18, 19, 20, 26, 27, 28, 29]
    for p_idx in target_indices:
        p = paras[p_idx]
        r = p[5] # child 5 has all gene names
        orig_text = ''.join(r.itertext())
        rPr = r.find('w:rPr', NS)
        
        pos = 0
        new_elements = []
        for m in pattern.finditer(orig_text):
            start, end = m.span()
            if start > pos:
                new_elements.append(make_plain_r(orig_text[pos:start], rPr))
            gene = m.group(1)
            new_elements.append(make_del(gene, rPr, cur_id))
            cur_id += 1
            new_elements.append(make_ins(gene, rPr, cur_id, italic=True))
            cur_id += 1
            pos = end
        if pos < len(orig_text):
            new_elements.append(make_plain_r(orig_text[pos:], rPr))
        
        # Replace r with new_elements
        r_idx = p.index(r)
        p.remove(r)
        for offset, el in enumerate(new_elements):
            p.insert(r_idx + offset, el)
        print(f"  P{p_idx}: Italicized gene symbols with tracked changes.")
    
    # 2. Punctuation fix in P36: 'repository.(' -> 'repository (' and ') Released' -> '). Released'
    p36 = paras[36]
    r1 = p36[1]
    t1 = r1.find('w:t', NS)
    if t1 is not None and t1.text.endswith("repository.("):
        # Change r1 text to end with "repository"
        t1.text = t1.text[:-2] # drop ".("
        rPr1 = r1.find('w:rPr', NS)
        r1_idx = p36.index(r1)
        # insert del of '.', plain ' ('
        del_dot = make_del(".", rPr1, cur_id)
        cur_id += 1
        plain_paren = make_plain_r(" (", rPr1)
        p36.insert(r1_idx + 1, del_dot)
        p36.insert(r1_idx + 2, plain_paren)
        print("  P36: Tracked deletion of period before '('")
    
    # Fix closing parenthesis in P36 child 5
    for c in p36:
        if c.tag.endswith('r'):
            t = c.find('w:t', NS)
            if t is not None and t.text and t.text.startswith(") Released"):
                rPr_c = c.find('w:rPr', NS)
                rem_text = t.text[1:] # " Released..."
                t.text = ")"
                c_idx = p36.index(c)
                ins_dot = make_ins(".", rPr_c, cur_id)
                cur_id += 1
                plain_rem = make_plain_r(rem_text, rPr_c)
                p36.insert(c_idx + 1, ins_dot)
                p36.insert(c_idx + 2, plain_rem)
                print("  P36: Tracked insertion of period after ')'")
                break
    
    tree.write(xml_path, encoding='utf-8', xml_declaration=True)
    repack_docx(temp_dir, docx_path)
    shutil.rmtree(temp_dir)
    print(f"Successfully updated Methods doc: {docx_path}\n")

# ==============================================================================
# Process Document 2: Mouse_flight_with_legends_SB.docx
# ==============================================================================
def process_sb_doc(docx_path):
    print(f"Processing SB doc: {docx_path}")
    temp_dir = tempfile.mkdtemp(prefix="sb_docx_")
    with zipfile.ZipFile(docx_path, 'r') as z:
        z.extractall(temp_dir)
    
    xml_path = os.path.join(temp_dir, 'word', 'document.xml')
    tree = etree.parse(xml_path)
    root = tree.getroot()
    body = root.find('w:body', NS)
    paras = body.findall('w:p', NS)
    
    ids = []
    for el in root.iter():
        for k, v in el.attrib.items():
            if k.endswith('id'):
                try:
                    ids.append(int(v))
                except ValueError:
                    pass
    cur_id = max(ids) + 1 if ids else 100
    
    # 1. P12: Mre11A -> Mre11a
    p12 = paras[12]
    c59 = p12[59]
    t59 = c59.find('w:t', NS)
    if t59 is not None and "Mre11A" in t59.text:
        rPr = c59.find('w:rPr', NS)
        idx = p12.index(c59)
        p12.remove(c59)
        del_el = make_del("Mre11A", rPr, cur_id)
        cur_id += 1
        ins_el = make_ins("Mre11a", rPr, cur_id, italic=True)
        cur_id += 1
        rem_r = make_plain_r(", Rad50, ", rPr)
        p12.insert(idx, del_el)
        p12.insert(idx + 1, ins_el)
        p12.insert(idx + 2, rem_r)
        print("  P12: Mre11A -> Mre11a applied.")
    
    # 2. P33: was -> were
    p33 = paras[33]
    c31 = p33[31]
    t31 = c31.find('w:t', NS)
    if t31 is not None and t31.text == "was":
        rPr = c31.find('w:rPr', NS)
        idx = p33.index(c31)
        p33.remove(c31)
        del_el = make_del("was", rPr, cur_id)
        cur_id += 1
        ins_el = make_ins("were", rPr, cur_id, italic=False)
        cur_id += 1
        p33.insert(idx, del_el)
        p33.insert(idx + 1, ins_el)
        print("  P33: was -> were applied.")
    
    # 3. P33: that increase -> that increases
    for child in list(p33):
        if child.tag.endswith('r'):
            t = child.find('w:t', NS)
            if t is not None and "that increase the robustness" in t.text:
                rPr = child.find('w:rPr', NS)
                idx = p33.index(child)
                p33.remove(child)
                # Split: ", that ", del "increase", ins "increases", " the robustness of NHEJ machinery and the "
                p33.insert(idx, make_plain_r(", that ", rPr))
                del_inc = make_del("increase", rPr, cur_id)
                cur_id += 1
                ins_inc = make_ins("increases", rPr, cur_id, italic=False)
                cur_id += 1
                post_r = make_plain_r(" the robustness of NHEJ machinery and the ", rPr)
                p33.insert(idx + 1, del_inc)
                p33.insert(idx + 2, ins_inc)
                p33.insert(idx + 3, post_r)
                print("  P33: that increase -> that increases applied.")
                break
    
    # 4. P33: ends ligation -> end ligation
    for child in list(p33):
        if child.tag.endswith('r'):
            t = child.find('w:t', NS)
            if t is not None and "and ends ligation" in t.text:
                rPr = child.find('w:rPr', NS)
                idx = p33.index(child)
                p33.remove(child)
                p33.insert(idx, make_plain_r(" and ", rPr))
                del_ends = make_del("ends", rPr, cur_id)
                cur_id += 1
                ins_end = make_ins("end", rPr, cur_id, italic=False)
                cur_id += 1
                post_lig = make_plain_r(" ligation", rPr)
                p33.insert(idx + 1, del_ends)
                p33.insert(idx + 2, ins_end)
                p33.insert(idx + 3, post_lig)
                print("  P33: ends ligation -> end ligation applied.")
                break
    
    # 5. P38: occured -> occurred
    p38 = paras[38]
    for child in list(p38):
        if child.tag.endswith('r'):
            t = child.find('w:t', NS)
            if t is not None and t.text == "occured":
                rPr = child.find('w:rPr', NS)
                idx = p38.index(child)
                p38.remove(child)
                del_occ = make_del("occured", rPr, cur_id)
                cur_id += 1
                ins_occ = make_ins("occurred", rPr, cur_id, italic=False)
                cur_id += 1
                p38.insert(idx, del_occ)
                p38.insert(idx + 1, ins_occ)
                print("  P38: occured -> occurred applied.")
                break
    
    # 6. P58: Ercc6 (plain) -> Ercc6 (italic)
    p58 = paras[58]
    for child in list(p58):
        if child.tag.endswith('r'):
            t = child.find('w:t', NS)
            if t is not None and "Ercc6 (" in t.text:
                rPr = child.find('w:rPr', NS)
                idx = p58.index(child)
                p58.remove(child)
                pre_text = t.text.split("Ercc6")[0] # " triangle), "
                post_text = t.text.split("Ercc6")[1] # " ("
                p58.insert(idx, make_plain_r(pre_text, rPr))
                del_ercc6 = make_del("Ercc6", rPr, cur_id)
                cur_id += 1
                ins_ercc6 = make_ins("Ercc6", rPr, cur_id, italic=True)
                cur_id += 1
                post_r = make_plain_r(post_text, rPr)
                p58.insert(idx + 1, del_ercc6)
                p58.insert(idx + 2, ins_ercc6)
                p58.insert(idx + 3, post_r)
                print("  P58: Ercc6 -> italic Ercc6 applied.")
                break
    
    tree.write(xml_path, encoding='utf-8', xml_declaration=True)
    repack_docx(temp_dir, docx_path)
    shutil.rmtree(temp_dir)
    print(f"Successfully updated SB doc: {docx_path}\n")

if __name__ == '__main__':
    methods_file = '/Users/gzy2520/Desktop/mouse_Methods.docx'
    sb_file = '/Users/gzy2520/Downloads/Mouse_flight_with_legends_SB.docx'
    process_methods_doc(methods_file)
    process_sb_doc(sb_file)
