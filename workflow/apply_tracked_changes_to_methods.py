#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Apply restrained, teacher-aligned tracked changes to Methods document XML:
1. Only cite Table S1, Table S2, Table S3, Table S4 where truly appropriate.
2. Use 'Table S' notation strictly, removing 'Supplementary'.
3. Remove redundant bullet citations and main figure table conversions.
"""

import copy
import lxml.etree as etree

tree = etree.parse('clean_unpack2/word/document.xml')
root = tree.getroot()
ns = {
    'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math'
}

cur_id = 350

def create_ins(text, author="gzy", date="2026-09-10T10:30:00Z"):
    global cur_id
    cur_id += 1
    ins = etree.Element(f"{{{ns['w']}}}ins", {
        f"{{{ns['w']}}}id": str(cur_id),
        f"{{{ns['w']}}}author": author,
        f"{{{ns['w']}}}date": date
    })
    r = etree.SubElement(ins, f"{{{ns['w']}}}r")
    rPr = etree.SubElement(r, f"{{{ns['w']}}}rPr")
    etree.SubElement(rPr, f"{{{ns['w']}}}rFonts", {
        f"{{{ns['w']}}}ascii": "Calibri",
        f"{{{ns['w']}}}hAnsi": "Calibri"
    })
    etree.SubElement(rPr, f"{{{ns['w']}}}color", {f"{{{ns['w']}}}val": "262626"})
    etree.SubElement(rPr, f"{{{ns['w']}}}sz", {f"{{{ns['w']}}}val": "20"})
    
    t = etree.SubElement(r, f"{{{ns['w']}}}t")
    if text.startswith(' ') or text.endswith(' '):
        t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    t.text = text
    return ins

def create_del(text, author="gzy", date="2026-09-10T10:30:00Z"):
    global cur_id
    cur_id += 1
    del_el = etree.Element(f"{{{ns['w']}}}del", {
        f"{{{ns['w']}}}id": str(cur_id),
        f"{{{ns['w']}}}author": author,
        f"{{{ns['w']}}}date": date
    })
    r = etree.SubElement(del_el, f"{{{ns['w']}}}r")
    rPr = etree.SubElement(r, f"{{{ns['w']}}}rPr")
    etree.SubElement(rPr, f"{{{ns['w']}}}rFonts", {
        f"{{{ns['w']}}}ascii": "Calibri",
        f"{{{ns['w']}}}hAnsi": "Calibri"
    })
    etree.SubElement(rPr, f"{{{ns['w']}}}color", {f"{{{ns['w']}}}val": "262626"})
    etree.SubElement(rPr, f"{{{ns['w']}}}sz", {f"{{{ns['w']}}}val": "20"})
    
    t = etree.SubElement(r, f"{{{ns['w']}}}delText")
    if text.startswith(' ') or text.endswith(' '):
        t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    t.text = text
    return del_el

paras = root.find('w:body', ns).findall('w:p', ns)

# 1. P2: Update ins id="2" from "(Supplementary Table S1)" to "(Table S1)"
p2 = paras[2]
for c in p2.findall('w:ins', ns):
    if c.attrib.get(f"{{{ns['w']}}}id") == "2":
        t = c.find('.//w:t', ns)
        if t is not None and "Supplementary Table S1" in t.text:
            idx = p2.index(c)
            p2.remove(c)
            p2.insert(idx, create_del(" (Supplementary Table S1)"))
            p2.insert(idx + 1, create_ins(" (Table S1)"))
            break

# 2. P3: Update ins id="8" from "(Supplementary Table S1)" to "(Table S1)"
p3 = paras[3]
for c in p3.findall('w:ins', ns):
    if c.attrib.get(f"{{{ns['w']}}}id") == "8":
        t = c.find('.//w:t', ns)
        if t is not None and "Supplementary Table S1" in t.text:
            idx = p3.index(c)
            p3.remove(c)
            p3.insert(idx, create_del(" (Supplementary Table S1)"))
            p3.insert(idx + 1, create_ins(" (Table S1)"))
            break

# 3. P4: Update "(Supplementary Table S1)" to "(Table S1)" in ins id="10"
p4 = paras[4]
for c in p4.findall('w:ins', ns):
    if c.attrib.get(f"{{{ns['w']}}}id") == "10":
        for t in c.findall('.//w:t', ns):
            if t.text and "Supplementary Table S1" in t.text:
                t.text = t.text.replace("Supplementary Table S1", "Table S1")
        break

# 4. P11: Delete gzy's ins id="73" "(Supplementary Table S5)" entirely via tracked deletion
p11 = paras[11]
for c in p11.findall('w:ins', ns):
    if c.attrib.get(f"{{{ns['w']}}}id") == "73":
        idx = p11.index(c)
        p11.remove(c)
        p11.insert(idx, create_del(" (Supplementary Table S5)"))
        break

# 5. P13: Add (Table S2) once at the end of the sentence before the period
p13 = paras[13]
last_r = p13.findall('w:r', ns)[-1]
t_last = last_r.find('w:t', ns)
if t_last is not None and t_last.text.endswith('.'):
    t_last.text = t_last.text[:-1]
    idx = p13.index(last_r)
    p13.insert(idx + 1, create_ins(" (Table S2)"))
    r_dot = etree.Element(f"{{{ns['w']}}}r")
    r_dot.append(copy.deepcopy(last_r.find('w:rPr', ns)))
    t_dot = etree.SubElement(r_dot, f"{{{ns['w']}}}t")
    t_dot.text = "."
    p13.insert(idx + 2, r_dot)

# 6. P14: Add (Table S2) after "11 in low-NES tissues"
p14 = paras[14]
for c in list(p14):
    if c.attrib.get(f"{{{ns['w']}}}id") == "137":
        idx = p14.index(c)
        p14.insert(idx + 1, create_ins(" (Table S2)"))
        break

# 7. P21: Add (Table S4) after "evaluated by one-way ANOVA across anatomical tissues."
p21 = paras[21]
for c in p21.findall('w:r', ns):
    t = c.find('w:t', ns)
    if t is not None and "across anatomical tissues." in (t.text or ""):
        t.text = t.text.replace("across anatomical tissues.", "across anatomical tissues")
        idx = p21.index(c)
        p21.insert(idx + 1, create_ins(" (Table S4)"))
        r_dot = etree.Element(f"{{{ns['w']}}}r")
        r_dot.append(copy.deepcopy(c.find('w:rPr', ns)))
        t_dot = etree.SubElement(r_dot, f"{{{ns['w']}}}t")
        t_dot.text = "."
        p21.insert(idx + 2, r_dot)
        break

# 8. P30: Add (Table S4) after "evaluated by one-way ANOVA."
p30 = paras[30]
for c in p30.findall('w:r', ns):
    t = c.find('w:t', ns)
    if t is not None and "evaluated by one-way ANOVA." in (t.text or ""):
        t.text = t.text.replace("evaluated by one-way ANOVA.", "evaluated by one-way ANOVA")
        idx = p30.index(c)
        p30.insert(idx + 1, create_ins(" (Table S4)"))
        r_dot = etree.Element(f"{{{ns['w']}}}r")
        r_dot.append(copy.deepcopy(c.find('w:rPr', ns)))
        t_dot = etree.SubElement(r_dot, f"{{{ns['w']}}}t")
        t_dot.text = "."
        p30.insert(idx + 2, r_dot)
        break

# 9. P33: Add (Table S3) after del id="230"
p33 = paras[33]
for c in list(p33):
    if c.attrib.get(f"{{{ns['w']}}}id") == "230":
        idx = p33.index(c)
        p33.insert(idx + 1, create_ins(" (Table S3)"))
        break

tree.write('clean_unpack2/word/document.xml', encoding='utf-8', xml_declaration=True)
print("Updated clean_unpack2/word/document.xml with restrained, teacher-aligned Table S1-S4 citations!")
