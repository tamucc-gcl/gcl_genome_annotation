"""Native BUSCO/OMArk summary parsers used by structural QC reporting."""
import re
from pathlib import Path

def parse_busco(path):
    text = Path(path).read_text(encoding='utf-8')
    result = {"busco_complete": 0, "busco_single": 0, "busco_duplicated": 0,
              "busco_fragmented": 0, "busco_missing": 0, "busco_total": 0}
    for key, pat in [("busco_complete","C:"), ("busco_single","S:"),
                     ("busco_duplicated","D:"), ("busco_fragmented","F:"),
                     ("busco_missing","M:")]:
        m = re.search(re.escape(pat) + r"([\d.]+)%", text)
        if m:
            result[key] = float(m.group(1))
    m = re.search(r"n:\s*(\d+)", text)
    if m:
        result["busco_total"] = int(m.group(1))
    return result

def parse_omark(path):
    text = Path(path).read_text(encoding='utf-8')
    result = {k: 0.0 for k in [
        "omark_single", "omark_duplicated", "omark_dup_unexpected",
        "omark_dup_expected", "omark_missing", "omark_consistent",
        "omark_inconsistent", "omark_contaminant", "omark_unknown",
        "omark_partial", "omark_fragmented", "omark_n_proteins",
    ]}

    for line in text.splitlines():
        ls = line.strip()

        # Completeness line: S:59.13%,D:24.67%[U:24.57%,E:0.10%],M:16.20%
        if ls.startswith("S:") and "D:" in ls and "M:" in ls and "%" in ls:
            m = re.search(r"S:([\d.]+)%", ls)
            if m: result["omark_single"] = float(m.group(1))
            m = re.search(r"D:([\d.]+)%", ls)
            if m: result["omark_duplicated"] = float(m.group(1))
            m = re.search(r"\[U:([\d.]+)%", ls)
            if m: result["omark_dup_unexpected"] = float(m.group(1))
            m = re.search(r",E:([\d.]+)%\]", ls)
            if m: result["omark_dup_expected"] = float(m.group(1))
            m = re.search(r",M:([\d.]+)%", ls)
            if m: result["omark_missing"] = float(m.group(1))

        # Consistency line: A:85.34%[P:15.56%,F:10.99%],I:6.10%,C:0.00%,U:8.56%
        if ls.startswith("A:") and "I:" in ls and "%" in ls:
            m = re.search(r"A:([\d.]+)%", ls)
            if m: result["omark_consistent"] = float(m.group(1))
            m = re.search(r"\[P:([\d.]+)%", ls)
            if m: result["omark_partial"] = float(m.group(1))
            m = re.search(r",F:([\d.]+)%\]", ls)
            if m: result["omark_fragmented"] = float(m.group(1))
            m = re.search(r",I:([\d.]+)%", ls)
            if m: result["omark_inconsistent"] = float(m.group(1))
            m = re.search(r",C:([\d.]+)%", ls)
            if m: result["omark_contaminant"] = float(m.group(1))
            m = re.search(r",U:([\d.]+)%", ls)
            if m: result["omark_unknown"] = float(m.group(1))

        # Protein count
        m = re.search(r"there are (\d+) proteins", ls)
        if m: result["omark_n_proteins"] = int(m.group(1))

    return result
