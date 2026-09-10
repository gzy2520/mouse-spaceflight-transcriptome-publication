#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Apply minimal tracked changes (<w:ins> and <w:del>) to user's Methods document XML,
strictly updating and adding Supplementary Table citations to match teacher's manuscript.
"""

import copy
import lxml.etree as etree

tree = etree.parse('clean_unpack/word/document.xml')
root = tree.getroot()
ns = {
    'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math'
}

cur_id = 300

def create_ins(text, author="gzy", date="2026-09-10T10:18:00Z", is_bold=False):
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
    if is_bold:
        etree.SubElement(rPr, f"{{{ns['w']}}}b")
    etree.SubElement(rPr, f"{{{ns['w']}}}color", {f"{{{ns['w']}}}val": "262626"})
    etree.SubElement(rPr, f"{{{ns['w']}}}sz", {f"{{{ns['w']}}}val": "20"})
    
    t = etree.SubElement(r, f"{{{ns['w']}}}t")
    if text.startswith(' ') or text.endswith(' '):
        t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    t.text = text
    return ins

def create_del(text, author="gzy", date="2026-09-10T10:18:00Z"):
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

# 1. P8: Add Table S5 after "NES matrix." and Table S6 after "complete cohort coverage."
p8 = paras[8]
for c in p8.findall('w:r', ns):
    t = c.find('w:t', ns)
    if t is not None and "mission-equal NES matrix." in (t.text or ""):
        orig_text = t.text
        s5_marker = "Figure 1c presents the complete 26 × 15 mission-equal NES matrix"
        s6_marker = "Supplementary Figure S1 displays the top 12 Hallmark pathways exhibiting the highest global positive and negative enrichment across missions with complete cohort coverage"
        
        idx_s5 = orig_text.find(s5_marker)
        idx_s6 = orig_text.find(s6_marker)
        
        before_s5 = orig_text[:idx_s5 + len(s5_marker)]
        between = orig_text[idx_s5 + len(s5_marker):idx_s6 + len(s6_marker)]
        after_s6 = orig_text[idx_s6 + len(s6_marker):]
        
        t.text = before_s5
        ins_s5 = create_ins(" (Supplementary Table S5)")
        
        r2 = etree.Element(f"{{{ns['w']}}}r")
        r2.append(copy.deepcopy(c.find('w:rPr', ns)))
        t2 = etree.SubElement(r2, f"{{{ns['w']}}}t")
        t2.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
        t2.text = between
        
        ins_s6 = create_ins(" (Supplementary Table S6)")
        
        r3 = etree.Element(f"{{{ns['w']}}}r")
        r3.append(copy.deepcopy(c.find('w:rPr', ns)))
        t3 = etree.SubElement(r3, f"{{{ns['w']}}}t")
        t3.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
        t3.text = after_s6
        
        idx = p8.index(c)
        p8.insert(idx + 1, ins_s5)
        p8.insert(idx + 2, r2)
        p8.insert(idx + 3, ins_s6)
        p8.insert(idx + 4, r3)
        break

# 2. P11: Replace gzy's ins id=73 "(Supplementary Table S5)" with del and ins of Table S8.
# And after "depicted in Figure 2", insert Table S7.
p11 = paras[11]
for c in p11.findall('w:ins', ns):
    if c.attrib.get(f"{{{ns['w']}}}id") == "73":
        idx = p11.index(c)
        p11.remove(c)
        p11.insert(idx, create_del(" (Supplementary Table S5)"))
        p11.insert(idx + 1, create_ins(" (Supplementary Table S8)"))
        break

for c in p11.findall('w:r', ns):
    t = c.find('w:t', ns)
    if t is not None and "depicted in Figure 2." in (t.text or ""):
        t.text = t.text.replace("depicted in Figure 2.", "depicted in Figure 2")
        idx = p11.index(c)
        ins_s7 = create_ins(" (Supplementary Table S7)")
        p11.insert(idx + 1, ins_s7)
        r_dot = etree.Element(f"{{{ns['w']}}}r")
        r_dot.append(copy.deepcopy(c.find('w:rPr', ns)))
        t_dot = etree.SubElement(r_dot, f"{{{ns['w']}}}t")
        t_dot.text = "."
        p11.insert(idx + 2, r_dot)
        break

# 3. P13: after del 113, 114, 115, 116 -> insert "; Supplementary Table S2"
p13 = paras[13]
ins_map_13 = {
    "113": "; Supplementary Table S2",
    "114": "; Supplementary Table S2",
    "115": "; Supplementary Table S2",
    "116": "; Supplementary Table S2"
}
for c in list(p13):
    del_id = c.attrib.get(f"{{{ns['w']}}}id")
    if del_id in ins_map_13:
        idx = p13.index(c)
        p13.insert(idx + 1, create_ins(ins_map_13[del_id]))

# 4. P14: after del 131, 133, 134, 137 -> insert Table S2
p14 = paras[14]
ins_map_14 = {
    "131": " (Supplementary Table S2),",
    "133": "(Supplementary Table S2)",
    "134": " (Supplementary Table S2)",
    "137": "(Supplementary Table S2), "
}
for c in list(p14):
    del_id = c.attrib.get(f"{{{ns['w']}}}id")
    if del_id in ins_map_14:
        idx = p14.index(c)
        p14.insert(idx + 1, create_ins(ins_map_14[del_id]))

# 5. P18, P19, P20: NHEJ, HR, A-EJ -> Table S4
for p_idx, del_target in [(18, "174"), (19, "177"), (20, "178")]:
    p = paras[p_idx]
    for c in list(p):
        if c.attrib.get(f"{{{ns['w']}}}id") == del_target:
            idx = p.index(c)
            p.insert(idx + 1, create_ins("Supplementary Table S4, ", is_bold=True))
            break

# 6. P21: ANOVA in DSB -> Table S4
p21 = paras[21]
for c in p21.findall('w:r', ns):
    t = c.find('w:t', ns)
    if t is not None and "across anatomical tissues." in (t.text or ""):
        t.text = t.text.replace("across anatomical tissues.", "across anatomical tissues")
        idx = p21.index(c)
        p21.insert(idx + 1, create_ins(" (Supplementary Table S4)"))
        r_dot = etree.Element(f"{{{ns['w']}}}r")
        r_dot.append(copy.deepcopy(c.find('w:rPr', ns)))
        t_dot = etree.SubElement(r_dot, f"{{{ns['w']}}}t")
        t_dot.text = "."
        p21.insert(idx + 2, r_dot)
        break

# 7. P22: Figure 4e -> Table S9
p22 = paras[22]
for c in list(p22):
    if c.attrib.get(f"{{{ns['w']}}}id") == "181":
        idx = p22.index(c)
        p22.insert(idx + 1, create_ins(" (Supplementary Table S9)"))
        break

# 8. P26, P27, P28, P29: BER, NER, MMR, FA -> Table S4
for p_idx, del_target in [(26, "213"), (27, "214"), (28, "217"), (29, "218")]:
    p = paras[p_idx]
    for c in list(p):
        if c.attrib.get(f"{{{ns['w']}}}id") == del_target:
            idx = p.index(c)
            p.insert(idx + 1, create_ins("Supplementary Table S4, ", is_bold=True))
            break

# 9. P30: ANOVA in SSB -> Table S4
p30 = paras[30]
for c in p30.findall('w:r', ns):
    t = c.find('w:t', ns)
    if t is not None and "evaluated by one-way ANOVA." in (t.text or ""):
        t.text = t.text.replace("evaluated by one-way ANOVA.", "evaluated by one-way ANOVA")
        idx = p30.index(c)
        p30.insert(idx + 1, create_ins(" (Supplementary Table S4)"))
        r_dot = etree.Element(f"{{{ns['w']}}}r")
        r_dot.append(copy.deepcopy(c.find('w:rPr', ns)))
        t_dot = etree.SubElement(r_dot, f"{{{ns['w']}}}t")
        t_dot.text = "."
        p30.insert(idx + 2, r_dot)
        break

# 10. P31: Figure 5f -> Table S10
p31 = paras[31]
for c in list(p31):
    if c.attrib.get(f"{{{ns['w']}}}id") == "221":
        idx = p31.index(c)
        p31.insert(idx + 1, create_ins(" (Supplementary Table S10)"))
        break

# 11. P33: Supplementary Figure S5 -> Table S3
p33 = paras[33]
for c in list(p33):
    if c.attrib.get(f"{{{ns['w']}}}id") == "230":
        idx = p33.index(c)
        p33.insert(idx + 1, create_ins("(Supplementary Table S3)"))
        break

tree.write('clean_unpack/word/document.xml', encoding='utf-8', xml_declaration=True)
print(f"Done! Clean document.xml updated with tracked changes up to ID {cur_id}.")
