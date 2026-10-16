"""Parse the coefficient tables of thesis appendix A (LaTeX source) into thesis_tables.json.

Usage: python3 parse_thesis.py /path/to/appendix_other_kernels.tex [out.json]

The tables write products like 5^63 and 3^714 meaning 5^6*3 and 3^7*14 (a single-digit exponent
followed by an integer factor); that convention is decoded in parse_expr().
Fractions are stored as strings "p/q".
"""
import json, re, sys
from fractions import Fraction as Fr

KEYS = {"Cubic spline": "M4", "Quartic spline": "M5", "Quintic spline": "M6",
        "Wendland $C^2$": "WC2", "Wendland $C^4$": "WC4", "Wendland $C^6$": "WC6"}

def parse_expr(s):
    """'5^63' -> 5^6*3, '2^3 3^8 7' -> 8*6561*7, '14' -> 14."""
    out = Fr(1)
    for tok in s.replace("\\pi", "").split():
        m = re.fullmatch(r"(\d+)\^(\d)(\d*)", tok)
        if m:
            out *= Fr(int(m[1])) ** int(m[2]) * (int(m[3]) if m[3] else 1)
        elif re.fullmatch(r"\d+", tok):
            out *= int(tok)
        else:
            raise ValueError("cannot parse %r" % s)
    return out

NUM = r"(-?)\s*(?:\\frac\{([^{}]*)\}\{([^{}]*)\}|(\d+))"
ENTRY = re.compile(r"\$?([ab])_(\{-?\d+\}|-?\d+)\s*=\s*" + NUM)
BREF = re.compile(r"a_0\s*=\s*B_(\d)")

def entries(text):
    out = {}
    for m in ENTRY.finditer(text):
        idx = int(m[2].strip("{}"))
        val = parse_expr(m[4]) / parse_expr(m[5]) if m[4] else Fr(int(m[6]))
        out[(m[1], idx)] = -val if m[3] else val
    return out

def body_rows(tab):
    lines = [l.split("%")[0] for l in tab.splitlines()]
    text = " ".join(lines).replace("\\hline", " ")
    text = text.split("\\label")[0] if False else text
    for chunk in text.split("\\\\"):
        chunk = chunk.strip()
        if chunk:
            yield [c.strip() for c in chunk.split("&")]

def main(src, dst):
    tex = open(src).read()
    out = {"source": "thesis appendix A (appendix_other_kernels.tex)", "A": {}, "2d": {}, "3d_a": {}, "3d_b": {}}
    for tab in re.findall(r"\\begin\{(?:sideways)?table\*?\}.*?\\end\{(?:sideways)?table\*?\}", tex, re.S):
        label = re.search(r"\\label\{table:([^}]*)\}", tab)[1]
        body = re.search(r"\\begin\{tabular\}\{[^}]*\}(.*?)\\end\{tabular\}", tab, re.S)[1]
        if label == "kernel-functions":
            for r in body_rows(body):
                if r[0] in KEYS:
                    a2 = re.search(r"\\frac\{([^{}]*)\}\{([^{}]*)\}", r[2]); a3 = re.search(r"\\frac\{([^{}]*)\}\{([^{}]*)\}", r[3])
                    f = lambda m: parse_expr(m[1]) / parse_expr(m[2])
                    out["A"][KEYS[r[0]]] = {"2": str(f(a2)), "3": str(f(a3))}
            continue
        kind = {"2d-IP-spline": "2d", "2d-IP-W": "2d", "3d-IP-spline": "3d_a", "3d-IP-W": "3d_a",
                "3d-Bn-spline": "3d_b", "3d-Bn-W": "3d_b"}[label]
        kern, region, qrange, prev = None, 0, 0, None
        for r in body_rows(body):
            if r[0] in KEYS:
                kern, region, qrange = KEYS[r[0]], 0, 0
            if kern is None or len(r) < 3 or "Range" in r[1]:
                continue
            cont = (r[0] == "" and r[1] == "" and (kind != "3d_b" or r[2] == ""))
            if kind == "3d_b":
                if r[1]:
                    qrange += 1
                if not cont:
                    k = int(re.search(r"B_(\d)", r[2])[1])
                    key = "%d,%d" % (qrange, k)
                    out["3d_b"].setdefault(kern, {})[key] = prev = {}
                for (t, idx), v in entries(" ".join(r)).items():
                    prev[str(idx)] = str(v)
            else:
                if not cont:
                    region += 1
                    d = out[kind].setdefault(kern, {})
                    prev = d[str(region)] = {}
                    if kind == "3d_a" and (m := BREF.search(" ".join(r))):
                        assert int(m[1]) == region, (kern, region, m[1])
                for (t, idx), v in entries(" ".join(r)).items():
                    prev[str(idx)] = str(v)
    json.dump(out, open(dst, "w"), indent=1)
    print("wrote", dst)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "thesis_tables.json")
