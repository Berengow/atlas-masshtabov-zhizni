# -*- coding: utf-8 -*-
"""Fetch structures, photos and silhouettes, then write build/data.json.
Run: python -I build/fetch.py  (from the project root)"""
import io, json, os, re, sys, time, html, urllib.request, urllib.parse, urllib.error
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from PIL import Image, ImageDraw, ImageFilter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "build"))
import content  # noqa

DL = os.path.join(ROOT, "downloads")
SITE = os.path.join(ROOT, "site")
for d in ("mol", "pdb", "rcsb", "wiki", "phylo", "meta"):
    os.makedirs(os.path.join(DL, d), exist_ok=True)
for d in ("img", "pdb"):
    os.makedirs(os.path.join(SITE, d), exist_ok=True)

UA = "ScaleAtlasBuilder/1.0 (educational, non-commercial)"
LOG = []

def log(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s, flush=True)

def get(url, binary=True, headers=None, tries=5):
    h = {"User-Agent": UA}
    if headers: h.update(headers)
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers=h)
            with urllib.request.urlopen(req, timeout=40) as r:
                data = r.read()
            return data if binary else data.decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code in (404, 400):
                return None
            ra = e.headers.get("Retry-After") if e.headers else None
            time.sleep(float(ra) if ra and ra.isdigit() else 4.0 * (t + 1))
        except Exception:
            time.sleep(1.5 * (t + 1))
    return None

def cached(path, fetch):
    if os.path.exists(path):
        with open(path, "rb") as f:
            d = f.read()
        return d if d != b"__NONE__" else None
    d = fetch()
    with open(path, "wb") as f:
        f.write(d if d else b"__NONE__")
    return d

# ───────────── molecules ─────────────
VDW = dict(H=1.2, C=1.7, N=1.55, O=1.52, S=1.8, P=1.8, Se=1.9, I=1.98, Fe=2.0, Co=2.0, Mg=1.73,
           Cl=1.75, Br=1.85, F=1.47, Na=2.27, K=2.75, Ca=2.31, Zn=1.39)
MASS = dict(H=1.008, C=12.011, N=14.007, O=15.999, S=32.06, P=30.974, Se=78.971, I=126.904, Fe=55.845,
            Co=58.933, Mg=24.305, Cl=35.45, Br=79.904, F=18.998, Na=22.990, K=39.098, Ca=40.078, Zn=65.38)

def parse_sdf(txt):
    lines = txt.splitlines()
    cnt = lines[3]
    na, nb = int(cnt[0:3]), int(cnt[3:6])
    atoms, bonds = [], []
    for l in lines[4:4 + na]:
        x, y, z = float(l[0:10]), float(l[10:20]), float(l[20:30])
        el = l[31:34].strip()
        el = el[0].upper() + el[1:].lower()
        atoms.append((el, x, y, z))
    for l in lines[4 + na:4 + na + nb]:
        a, b, o = int(l[0:3]), int(l[3:6]), int(l[6:9])
        bonds.append((a - 1, b - 1, o))
    return atoms, bonds

def mol_fetch(it):
    path = os.path.join(DL, "mol", it["id"] + ".sdf")
    def f():
        if it.get("pc"):
            q = urllib.parse.quote(it["pc"])
            d = get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/SDF?record_type=3d")
            time.sleep(0.25)
            if d: return d
        if it.get("lig"):
            for kind in ("ideal", "model"):
                d = get(f"https://files.rcsb.org/ligands/download/{it['lig']}_{kind}.sdf")
                if d and b"V2000" in d[:2000]:
                    # some ideal files have all-zero coords; check
                    a, _ = parse_sdf(d.decode())
                    if any(abs(x) + abs(y) + abs(z) > 0.01 for _, x, y, z in a):
                        return d
        if it.get("pc"):
            q = urllib.parse.quote(it["pc"])
            d = get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/SDF")
            time.sleep(0.25)
            if d: return d
        return None
    d = cached(path, f)
    if not d:
        log("MOL MISSING", it["id"])
        return None
    txt = d.decode("utf-8", "replace")
    atoms, bonds = parse_sdf(txt)
    els = [a[0] for a in atoms]
    X = np.array([[a[1], a[2], a[3]] for a in atoms])
    X = X - X.mean(0)
    if len(X) >= 3:
        cov = np.cov(X.T)
        w, V = np.linalg.eigh(cov)
        V = V[:, ::-1]
        if np.linalg.det(V) < 0: V[:, 2] *= -1
        X = X @ V
    r = np.array([VDW.get(e, 1.8) for e in els])
    xmin, xmax = (X[:, 0] - r).min(), (X[:, 0] + r).max()
    ymin, ymax = (X[:, 1] - r).min(), (X[:, 1] + r).max()
    X[:, 0] -= (xmin + xmax) / 2
    X[:, 1] -= (ymin + ymax) / 2
    wA, hA = xmax - xmin, ymax - ymin
    # formula (Hill)
    from collections import Counter
    c = Counter(els)
    order = (["C", "H"] if "C" in c else []) + sorted(k for k in c if k not in (("C", "H") if "C" in c else ()))
    sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    formula = "".join(k + (str(c[k]).translate(sub) if c[k] > 1 else "") for k in order)
    mass = sum(MASS.get(e, 0) for e in els)
    two_d = bool(np.allclose(X[:, 2], 0))
    return dict(e=els, x=[round(float(v), 2) for v in X.flatten()], b=[v for bb in bonds for v in bb],
                w=float(wA), h=float(hA), formula=formula, mass=round(mass, 2), flat=two_d)

# ───────────── PDB ─────────────
NUC = {"DA", "DT", "DG", "DC", "A", "U", "G", "C", "DU"}
def pdb_fetch(pid):
    raw = cached(os.path.join(DL, "pdb", pid + ".pdb"),
                 lambda: get(f"https://files.rcsb.org/download/{pid}.pdb"))
    if not raw:
        log("PDB MISSING", pid); return None
    lines = raw.decode("utf-8", "replace").splitlines()
    keep = []
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            res = l[17:20].strip()
            if res in ("HOH", "WAT", "DOD"): continue
            el = l[76:78].strip() or l[12:16].strip()[0]
            if el in ("H", "D"): continue
            keep.append(l)
        elif l.startswith(("TER", "ENDMDL")):
            keep.append(l)
            if l.startswith("ENDMDL"): break
    txt = "\n".join(keep)
    if len(txt) > 900_000:
        bb = []
        for l in keep:
            if l.startswith("ATOM"):
                name = l[12:16].strip(); res = l[17:20].strip()
                if res in NUC or name in ("N", "CA", "C", "O"):
                    bb.append(l)
            elif l.startswith("TER"):
                bb.append(l)
            elif l.startswith("HETATM"):
                bb.append(l)
        txt = "\n".join(bb)
    txt += "\nEND\n"
    with open(os.path.join(SITE, "pdb", pid + ".pdb"), "w", encoding="utf-8") as f:
        f.write(txt)
    return "pdb/" + pid + ".pdb", len(txt)

def remove_white(im):
    im = im.convert("RGB")
    W, H = im.size
    a = np.asarray(im).astype(np.int16)
    near = (a.min(2) > 228)
    # flood fill from border over near-white pixels
    mask = Image.fromarray((near * 255).astype(np.uint8))
    seed_img = Image.new("L", (W + 2, H + 2), 255)
    seed_img.paste(mask, (1, 1))
    ImageDraw.floodfill(seed_img, (0, 0), 128, thresh=0)  # fills the near-white background connected to the border
    bg = (np.asarray(seed_img)[1:-1, 1:-1] == 128)
    mn = a.min(2).astype(np.float32)
    soft = np.clip((246.0 - mn) / (246.0 - 212.0), 0, 1) * 255  # interior white gaps fade out too
    alpha = np.where(bg, 0, soft).astype(np.uint8)
    am = Image.fromarray(alpha).filter(ImageFilter.GaussianBlur(0.9))
    rgba = im.convert("RGBA"); rgba.putalpha(am)
    bbox = am.point(lambda v: 255 if v > 20 else 0).getbbox()
    if bbox: rgba = rgba.crop(bbox)
    return rgba

def rcsb_image(it, pid):
    p = pid.lower()
    def f():
        for suf in ("assembly-1", "model-1"):
            d = get(f"https://cdn.rcsb.org/images/structures/{p}_{suf}.jpeg")
            if d: return d
        return None
    d = cached(os.path.join(DL, "rcsb", p + ".jpeg"), f)
    if not d:
        log("RCSB IMG MISSING", pid); return None
    im = remove_white(Image.open(io.BytesIO(d)))
    im.thumbnail((720, 720), Image.LANCZOS)
    out = f"img/p_{it['id']}.webp"
    im.save(os.path.join(SITE, out), "WEBP", quality=86, method=6)
    return dict(src=out, w=im.size[0], h=im.size[1])

# ───────────── Wikimedia photos ─────────────
def strip_html(s):
    s = re.sub(r"<[^>]+>", "", s or "")
    return html.unescape(re.sub(r"\s+", " ", s)).strip()

def commons_info(fname):
    q = urllib.parse.urlencode(dict(action="query", titles="File:" + fname, prop="imageinfo",
                                    iiprop="url|extmetadata|size|mime", iiurlwidth=720, format="json"))
    t = get("https://commons.wikimedia.org/w/api.php?" + q, binary=False)
    if not t: return None
    pages = json.loads(t)["query"]["pages"]
    pg = next(iter(pages.values()))
    if "imageinfo" not in pg: return None
    ii = pg["imageinfo"][0]
    md = ii.get("extmetadata", {})
    return dict(file=fname, thumb=ii.get("thumburl") or ii["url"], page=ii.get("descriptionurl"),
                mime=ii.get("mime"),
                artist=strip_html(md.get("Artist", {}).get("value", ""))[:140],
                license=strip_html(md.get("LicenseShortName", {}).get("value", "")),
                desc=strip_html(md.get("ImageDescription", {}).get("value", ""))[:200])

def wiki_lead(title):
    q = urllib.parse.urlencode(dict(action="query", titles=title, prop="pageimages", piprop="name",
                                    pilicense="free", format="json", redirects=1))
    t = get("https://en.wikipedia.org/w/api.php?" + q, binary=False)
    if not t: return None
    pages = json.loads(t)["query"]["pages"]
    pg = next(iter(pages.values()))
    return pg.get("pageimage")

def commons_search(query):
    q = urllib.parse.urlencode(dict(action="query", list="search", srnamespace=6, srsearch=query,
                                    srlimit=8, format="json"))
    t = get("https://commons.wikimedia.org/w/api.php?" + q, binary=False)
    if not t: return None
    for r in json.loads(t)["query"]["search"]:
        name = r["title"].split(":", 1)[1]
        if re.search(r"\.(jpe?g|png|tif?f)$", name, re.I):
            return name
    return None

def photo(it):
    meta_p = os.path.join(DL, "meta", "w_" + it["id"] + ".json")
    if os.path.exists(meta_p):
        meta = json.load(open(meta_p, encoding="utf-8"))
    else:
        fname = wiki_lead(it["w"]) if it.get("w") else None
        if not fname and it.get("cs"):
            fname = commons_search(it["cs"])
        meta = commons_info(fname) if fname else None
        if meta and it.get("w"): meta["wiki"] = "https://en.wikipedia.org/wiki/" + urllib.parse.quote(it["w"].replace(" ", "_"))
        json.dump(meta, open(meta_p, "w", encoding="utf-8"), ensure_ascii=False)
    if not meta:
        log("PHOTO MISSING", it["id"]); return None
    raw = cached(os.path.join(DL, "wiki", it["id"] + ".bin"), lambda: get(meta["thumb"]))
    if not raw:
        log("PHOTO DL FAIL", it["id"]); return None
    try:
        im = Image.open(io.BytesIO(raw))
        if im.mode in ("RGBA", "LA", "P"):
            im = im.convert("RGBA"); bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
            bg.alpha_composite(im); im = bg
        im = im.convert("RGB")
        im.thumbnail((720, 720), Image.LANCZOS)
        out = f"img/w_{it['id']}.jpg"
        im.save(os.path.join(SITE, out), "JPEG", quality=80, optimize=True, progressive=True)
    except Exception as e:
        log("PHOTO BAD", it["id"], e); return None
    return dict(src=out, w=im.size[0], h=im.size[1], artist=meta.get("artist"), license=meta.get("license"),
                page=meta.get("page"), wiki=meta.get("wiki"), file=meta.get("file"))

# ───────────── PhyloPic ─────────────
BUILD = None
def phylo(it):
    global BUILD
    meta_p = os.path.join(DL, "meta", "pp_" + it["id"] + ".json")
    if os.path.exists(meta_p):
        meta = json.load(open(meta_p, encoding="utf-8"))
    else:
        meta = None
        for name in it["pp"]:
            q = urllib.parse.urlencode(dict(build=BUILD, filter_name=name.lower(), embed_items="true",
                                            embed_primaryImage="true", page=0))
            t = get("https://api.phylopic.org/nodes?" + q, binary=False,
                    headers={"Accept": "application/vnd.phylopic.v2+json"})
            if not t: continue
            d = json.loads(t)
            items = d.get("_embedded", {}).get("items", [])
            for n in items:
                img = n.get("_embedded", {}).get("primaryImage")
                if not img: continue
                L = img["_links"]
                rasters = sorted(L.get("rasterFiles", []), key=lambda r: int(r["sizes"].split("x")[0]))
                pick = next((r for r in rasters if int(r["sizes"].split("x")[0]) >= 1000), rasters[-1] if rasters else None)
                if not pick: continue
                meta = dict(src=pick["href"], name=name, node=n.get("names", [[{}]])[0][0].get("text") if n.get("names") else name,
                            license=L.get("license", {}).get("href"), contributor=L.get("contributor", {}).get("title"),
                            attribution=img.get("attribution"), page="https://www.phylopic.org/images/" + img["uuid"] if img.get("uuid") else None)
                break
            if meta: break
        json.dump(meta, open(meta_p, "w", encoding="utf-8"), ensure_ascii=False)
    if not meta:
        log("PHYLO MISSING", it["id"]); return None
    raw = cached(os.path.join(DL, "phylo", it["id"] + ".png"), lambda: get(meta["src"]))
    if not raw: return None
    im = Image.open(io.BytesIO(raw)).convert("RGBA")
    a = im.split()[3]
    bbox = a.getbbox()
    if bbox: im = im.crop(bbox)
    im.thumbnail((1000, 1000), Image.LANCZOS)
    # store as white silhouette with alpha (tinted at runtime)
    a = im.split()[3]
    white = Image.new("RGBA", im.size, (255, 255, 255, 0)); white.putalpha(a)
    out = f"img/s_{it['id']}.webp"
    white.save(os.path.join(SITE, out), "WEBP", lossless=True, quality=100, method=6)
    lic = meta.get("license") or ""
    lic_s = "CC0" if "zero" in lic else ("Public Domain" if "publicdomain/mark" in lic else
             ("CC BY" + (("-NC" if "-nc" in lic else "") + ("-SA" if "-sa" in lic else ""))) if "by" in lic else lic)
    return dict(src=out, w=im.size[0], h=im.size[1], who=meta.get("attribution") or meta.get("contributor"),
                license=lic_s, page=meta.get("page"))

def main():
    global BUILD
    t = get("https://api.phylopic.org/", binary=False, headers={"Accept": "application/vnd.phylopic.v2+json"})
    BUILD = json.loads(t)["build"] if t else 562
    items = content.ITEMS
    out = {}

    def work(it):
        r = {}
        try:
            if it.get("pc") or it.get("lig"):
                r["mol"] = mol_fetch(it)
            if it.get("pdb") and not it.get("pdbimg"):
                pr = pdb_fetch(it["pdb"])
                if pr: r["pdbfile"] = pr[0]
            pid = it.get("pdb") or it.get("img3d")
            if pid:
                r["pimg"] = rcsb_image(it, pid)
            if it.get("w") or it.get("cs"):
                r["photo"] = photo(it)
            if it.get("pp"):
                r["sil"] = phylo(it)
        except Exception as e:
            import traceback; traceback.print_exc()
            log("ERR", it["id"], e)
        return it["id"], r

    mol_items = [i for i in items if i.get("pc") or i.get("lig")]
    other = [i for i in items if not (i.get("pc") or i.get("lig"))]
    with ThreadPoolExecutor(3) as ex:
        for k, r in ex.map(work, mol_items): out[k] = r
    with ThreadPoolExecutor(2) as ex:
        for k, r in ex.map(work, other): out[k] = r
    json.dump(out, open(os.path.join(ROOT, "build", "fetched.json"), "w", encoding="utf-8"), ensure_ascii=False)
    print("DONE", len(out))
    print("\n".join(l for l in LOG if "MISSING" in l or "ERR" in l or "BAD" in l or "FAIL" in l))

if __name__ == "__main__":
    main()
