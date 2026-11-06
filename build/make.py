# -*- coding: utf-8 -*-
"""Assemble site/index.html from content.py + fetched assets.
Run from the project root: python -I build/make.py [--refresh]"""
import json, os, re, sys, io

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "build"))
import content, fetch  # noqa

SMALL = {"h2o", "o2", "n2", "co2", "no", "h2s"}
AR = dict(atom=1, icosa=1, env=1, giantv=1, vesicle=1, pore=1, lyso=1, perox=1, myco=1.3, square=1.05, wbc=1, mega=1,
          adipo=1, egg=1, cocco=1, chlorella=1, diatom=1, radiolaria=1, noctiluca=1, volvox=1, platelet=1.3, tree=.25,
          rbc=1.7, phage=.45)

def main(refresh):
    fetched = json.load(open(os.path.join(ROOT, "build", "fetched.json"), encoding="utf-8"))
    byid = {i["id"]: i for i in content.ITEMS}
    if refresh:
        fetch.BUILD = 562
        for it in content.ITEMS:
            r = fetched.setdefault(it["id"], {})
            pid = it.get("pdb") or it.get("img3d")
            if pid:
                r["pimg"] = fetch.rcsb_image(it, pid)
            if it["id"] in SMALL:
                r["mol"] = fetch.mol_fetch(it)
            if it["id"] == "ebola":
                r["photo"] = fetch.photo({k: v for k, v in it.items() if k != "w"})
        # crisp vector silhouette where PhyloPic only had a tiny raster
        for iid, r in fetched.items():
            s = r.get("sil")
            if s and max(s["w"], s["h"]) < 400:
                meta = json.load(open(os.path.join(fetch.DL, "meta", "pp_" + iid + ".json"), encoding="utf-8"))
                vurl = re.sub(r"/raster/.*$", "/vector.svg", meta["src"])
                svg = fetch.get(vurl)
                if svg:
                    txt = svg.decode("utf-8")
                    m = re.search(r'viewBox="([\d.\s-]+)"', txt)
                    vb = [float(v) for v in m.group(1).split()] if m else [0, 0, s["w"], s["h"]]
                    txt = re.sub(r"<svg ", '<svg width="%d" height="%d" ' % (round(vb[2]), round(vb[3])), txt, count=1) if 'width="' not in txt[:400] else txt
                    txt = re.sub(r"<\?xml[^>]*\?>\s*", "", txt); txt = re.sub(r"<!DOCTYPE[^>]*>\s*", "", txt, flags=re.I)
                    out = "img/s_%s.svg" % iid
                    open(os.path.join(ROOT, "site", out), "w", encoding="utf-8").write(txt)
                    s.update(src=out, w=vb[2], h=vb[3])
        json.dump(fetched, open(os.path.join(ROOT, "build", "fetched.json"), "w", encoding="utf-8"), ensure_ascii=False)

    items = []
    for it in content.ITEMS:
        f = fetched.get(it["id"], {})
        o = {k: it[k] for k in ("id", "n", "en", "c", "d", "f") if k in it}
        if it.get("sd"): o["sd"] = it["sd"]
        if it.get("p"): o["p"] = it["p"]
        if it.get("pc"): o["pc"] = it["pc"]
        if it.get("lig"): o["lig"] = it["lig"]
        if it.get("w"): o["wiki"] = "https://en.wikipedia.org/wiki/" + it["w"].replace(" ", "_")
        v = dict(it.get("v") or {})
        s = it.get("s")
        mol = f.get("mol")
        pimg = f.get("pimg")
        sil = f.get("sil")
        ph = f.get("photo")
        if mol:
            if it["id"] in SMALL: mol["flat"] = False
            o["m"] = mol; v = {"t": "mol"}
            r2 = lambda x: float("%.2g" % x)
            w, h = r2(mol["w"]) * 1e-10, r2(mol["h"]) * 1e-10
        elif pimg and (not it.get("v") or it.get("pdbimg")):
            o["pi"] = pimg; v = {"t": "img"}
            ar = pimg["w"] / pimg["h"]
            w, h = (s, s / ar) if ar >= 1 else (s * ar, s)
        elif sil:
            o["sil"] = sil; v = {"t": "sil"}
            ar = sil["w"] / sil["h"]
            if it.get("dim") == "h": w, h = s * ar, s
            else: w, h = (s, s / ar) if ar >= 1 else (s * ar, s)
        else:
            ar = it.get("ar") or AR.get(v.get("t"), 1)
            if v.get("t") == "coccus" and v.get("arr") == "cluster": ar = 1.1
            if v.get("t") == "pollen" and v.get("kind") == "smooth": ar = 1.6
            w, h = (s, s / ar) if ar >= 1 else (s * ar, s)
        if pimg and "pi" not in o and it.get("pdbimg"): o["pi"] = pimg
        o["v"] = v; o["w"] = w; o["h"] = h
        pid = it.get("pdb") or it.get("img3d")
        if pid: o["pdbid"] = pid
        if f.get("pdbfile"):
            o["pdb"] = f["pdbfile"]
            txt = open(os.path.join(ROOT, "site", f["pdbfile"]), encoding="utf-8").read(400000)
            o["pdbbb"] = " CB " not in txt
        # photos: proteins with a 3D model don't need an extra diagram
        if ph and not (it["c"] == "prot" and f.get("pdbfile")):
            ph = dict(ph)
            art = re.split(r"\.\s+(?:When|If|Please|Use|Attribution|Credit)", ph.get("artist") or "")[0]
            ph["artist"] = (art[:77] + "…") if len(art) > 80 else art
            o["ph"] = ph
        items.append(o)

    order = {c[0]: i for i, c in enumerate(content.CHAPTERS)}
    items.sort(key=lambda o: (order[o["c"]], max(o["w"], o["h"])))
    chapters = [dict(id=c[0], name=c[1], blurb=c[2]) for c in content.CHAPTERS]
    data = json.dumps(dict(chapters=chapters, items=items), ensure_ascii=False, separators=(",", ":"))
    data = data.replace("</", "<\\/")
    tpl = open(os.path.join(ROOT, "build", "template.html"), encoding="utf-8").read()
    html = tpl.replace("__DATA__", data)
    open(os.path.join(ROOT, "site", "index.html"), "w", encoding="utf-8").write(html)
    write_docs(html, items)
    print("items", len(items), "html KB", len(html.encode()) // 1024)
    print("no photo:", [o["id"] for o in items if "ph" not in o and o["v"]["t"] not in ("mol", "atom") and "pdb" not in o])

def write_docs(html, items):
    """Standalone build for GitHub Pages in docs/: full HTML document plus only the assets it references."""
    import shutil
    docs = os.path.join(ROOT, "docs")
    for sub in ("img", "pdb"):
        os.makedirs(os.path.join(docs, sub), exist_ok=True)
    head, body = html.split("</style>", 1)
    page = ('<!doctype html>\n<html lang="ru">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
            '<meta property="og:title" content="Атлас масштабов жизни">\n'
            '<meta property="og:description" content="От атома водорода до гигантской грибницы: 262 объекта в одном масштабе, с 3D-моделями и фотографиями.">\n'
            '<style>:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}'
            'img{max-width:100%}[hidden]{display:none!important}</style>\n'
            + head + "</style>\n</head>\n<body>" + body + "\n</body>\n</html>\n")
    open(os.path.join(docs, "index.html"), "w", encoding="utf-8").write(page)
    open(os.path.join(docs, ".nojekyll"), "w").close()
    used = set(re.findall(r'"((?:img|pdb)/[^"]+)"', html))
    for rel in sorted(used):
        src = os.path.join(ROOT, "site", rel)
        if os.path.exists(src):  # fresh download wins; otherwise keep the committed copy in docs/
            shutil.copy2(src, os.path.join(docs, rel))
    for sub in ("img", "pdb"):
        for name in os.listdir(os.path.join(docs, sub)):
            if sub + "/" + name not in used:
                os.remove(os.path.join(docs, sub, name))
    # credits for every third-party image
    L = ["# Источники изображений и данных", "",
         "Фотографии — Wikimedia Commons, силуэты — PhyloPic, изображения белков и вирусов — RCSB PDB (rcsb.org), "
         "3D-координаты молекул — PubChem и RCSB PDB Chemical Component Dictionary. "
         "Лицензия каждого изображения указана ниже; при повторном использовании соблюдайте её условия.", "",
         "## Фотографии (Wikimedia Commons)", "", "| Объект | Файл | Автор | Лицензия |", "|---|---|---|---|"]
    cell = lambda v: (v or "—").replace("|", "/").replace("\n", " ")
    for o in items:
        ph = o.get("ph")
        if ph:
            L.append("| %s | [%s](%s) | %s | %s |" % (o["n"], cell(ph.get("file")), ph.get("page") or "", cell(ph.get("artist")), cell(ph.get("license"))))
    L += ["", "## Силуэты (PhyloPic)", "", "| Объект | Автор | Лицензия | Страница |", "|---|---|---|---|"]
    for o in items:
        s = o.get("sil")
        if s:
            L.append("| %s | %s | %s | %s |" % (o["n"], cell(s.get("who")), cell(s.get("license")), s.get("page") or "—"))
    L += ["", "## Структуры (RCSB PDB)", "", "| Объект | PDB ID |", "|---|---|"]
    for o in items:
        if o.get("pdbid"):
            L.append("| %s | [%s](https://www.rcsb.org/structure/%s) |" % (o["n"], o["pdbid"], o["pdbid"]))
    open(os.path.join(ROOT, "CREDITS.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("docs assets", len(used))

if __name__ == "__main__":
    main("--refresh" in sys.argv)
