# -*- coding: utf-8 -*-
"""Export deck.html to PDF and PPTX.

    python3 export.py          # both
    python3 export.py pdf
    python3 export.py pptx

The PDF is printed straight from Chrome, so its text stays vector and
selectable. The PPTX is one full-bleed 2560x1440 image per slide: PowerPoint
cannot represent this deck's inline SVG and CSS layout natively, so anything
else would be a lossy re-drawing rather than the deck itself.

Neither output needs a third-party library - python-pptx is not a dependency,
the OOXML is written directly. Both outputs are regenerated artifacts and are
gitignored; deck.html is the source.
"""
import os, re, sys, glob, shutil, zipfile, datetime, subprocess

HERE   = os.path.dirname(os.path.abspath(__file__))
DECK   = os.path.join(HERE, "deck.html")
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
TMP    = os.path.join(HERE, ".export_tmp")
W, H   = 12192000, 6858000          # 13.333in x 7.5in in EMU - standard 16:9

HEAD = ('<!doctype html><html data-theme="light"><head><meta charset=utf8>'
        '<style>:root{color-scheme:light}body{margin:0;padding:0;'
        'font:14px -apple-system,BlinkMacSystemFont,sans-serif;background:#faf9f5;color:#141413}'
        'img{max-width:100%}</style></head><body>')

# Printing needs three things the on-screen deck does not do: every slide laid
# out as its own page-sized block, backgrounds actually reproduced, and the
# entrance animation suppressed so a half-played frame is never captured. The
# clone drops the #bar/#count ids, so those id-keyed rules are restated here.
PRINT = '''
<style>
  @page { size: 13.333in 7.5in; margin: 0; }
  html, body { margin:0 !important; padding:0 !important; background:#fff !important; }
  *, *::before, *::after { -webkit-print-color-adjust: exact !important;
                           print-color-adjust: exact !important; }
  .pstage { width:1280px; height:720px; position:relative; overflow:hidden;
            background:var(--paper); break-after:page; page-break-after:always; }
  .pstage:last-child { break-after:auto; page-break-after:auto; }
  .slide.is-on > * { animation:none !important; }
  #fit { display:none !important; }
  .pbar { flex:1; display:flex; gap:3px; }
  .pbar button { flex:1; height:3px; border:0; padding:0; background:var(--rule);
                 border-radius:2px; -webkit-appearance:none; appearance:none; }
  .pbar button.done { background:var(--muted); }
  .pbar button.cur  { background:var(--accent); height:5px; margin-top:-1px; }
  .pcount { font-size:12px; color:var(--muted); min-width:70px; text-align:right; }
</style>
<script>
window.addEventListener('load', function(){
  var stage = document.getElementById('stage');
  var n = document.querySelectorAll('.slide').length;
  var holder = document.createElement('div');
  for (var i = 0; i < n; i++) {
    var c = stage.cloneNode(true);
    c.removeAttribute('id'); c.className = 'pstage'; c.style.transform = 'none';
    var sl = c.querySelectorAll('.slide');
    for (var k = 0; k < sl.length; k++) sl[k].classList.toggle('is-on', k === i);
    var act = sl[i].getAttribute('data-act');
    c.querySelectorAll('.act').forEach(function(a){
      a.classList.toggle('on', a.getAttribute('data-act') === act); });
    var cnt = c.querySelector('#count');
    if (cnt) { cnt.removeAttribute('id'); cnt.className = 'mn pcount';
               cnt.textContent = (i + 1) + ' / ' + n; }
    var bar = c.querySelector('#bar');
    if (bar) { bar.removeAttribute('id'); bar.className = 'pbar';
      [].forEach.call(bar.children, function(b, k){
        b.classList.toggle('cur', k === i); b.classList.toggle('done', k < i); }); }
    var hint = c.querySelector('#hint'); if (hint) hint.remove();
    holder.appendChild(c);
  }
  document.body.appendChild(holder);
});
</script>'''


def deck_html():
    if not os.path.exists(DECK):
        sys.exit("deck.html not found - run build.py first")
    return open(DECK, encoding="utf-8").read()


def chrome(*args):
    subprocess.run([CHROME, "--headless", "--disable-gpu", *args], capture_output=True)


def export_pdf():
    out = os.path.join(HERE, "PPIL4_deck.pdf")
    os.makedirs(TMP, exist_ok=True)
    p = os.path.join(TMP, "print.html")
    open(p, "w", encoding="utf-8").write(HEAD + deck_html() + PRINT + "</body></html>")
    chrome("--no-pdf-header-footer", "--print-to-pdf=" + out,
           "--virtual-time-budget=15000", "file://" + p)
    print("  PDF   %s  (%.1f MB)" % (out, os.path.getsize(out) / 1048576))
    return out


def render_slides():
    """One 2560x1440 PNG per slide, light theme, entrance animation suppressed."""
    html = deck_html()
    n = len(re.findall(r'<section class="slide"', html))
    d = os.path.join(TMP, "img")
    os.makedirs(d, exist_ok=True)
    head = HEAD.replace("img{max-width:100%}",
                        "img{max-width:100%}.slide.is-on > *{animation:none !important}")
    for k in range(n):
        p = os.path.join(d, "f%02d.html" % k)
        open(p, "w", encoding="utf-8").write(
            head + html.replace("  go(0);", "  go(%d);" % k) + "</body></html>")
        chrome("--hide-scrollbars", "--force-device-scale-factor=2",
               "--screenshot=" + os.path.join(d, "slide%02d.png" % (k + 1)),
               "--window-size=1280,720", "--virtual-time-budget=6000", "file://" + p)
        os.remove(p)
    return sorted(glob.glob(os.path.join(d, "slide*.png")))


NS = ('xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
      'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
      'xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"')
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/"
XML = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
EMPTY_TREE = ('<p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/>'
              '</p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/>'
              '<a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr></p:spTree>')


def _rels(pairs):
    return (XML + '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            + "".join('<Relationship Id="%s" Type="%s" Target="%s"/>' % p for p in pairs)
            + '</Relationships>')


def _theme():
    c = lambda n, v: '<a:%s><a:srgbClr val="%s"/></a:%s>' % (n, v, n)
    scheme = (c("dk1", "141A19") + c("lt1", "F5F7F6") + c("dk2", "3A4442") + c("lt2", "EDF0EF")
              + c("accent1", "12695E") + c("accent2", "0D8A76") + c("accent3", "4A3AA7")
              + c("accent4", "9C6B15") + c("accent5", "A8365C") + c("accent6", "626D6B")
              + c("hlink", "12695E") + c("folHlink", "626D6B"))
    font = ('<a:majorFont><a:latin typeface="Helvetica Neue"/><a:ea typeface=""/><a:cs typeface=""/></a:majorFont>'
            '<a:minorFont><a:latin typeface="Helvetica Neue"/><a:ea typeface=""/><a:cs typeface=""/></a:minorFont>')
    fill = '<a:solidFill><a:schemeClr val="phClr"/></a:solidFill>'
    ln = ('<a:ln w="9525" cap="flat" cmpd="sng" algn="ctr"><a:solidFill>'
          '<a:schemeClr val="phClr"/></a:solidFill><a:prstDash val="solid"/></a:ln>')
    fmt = ('<a:fmtScheme name="Office"><a:fillStyleLst>' + fill * 3 + '</a:fillStyleLst>'
           '<a:lnStyleLst>' + ln * 3 + '</a:lnStyleLst><a:effectStyleLst>'
           + '<a:effectStyle><a:effectLst/></a:effectStyle>' * 3
           + '</a:effectStyleLst><a:bgFillStyleLst>' + fill * 3 + '</a:bgFillStyleLst></a:fmtScheme>')
    return (XML + '<a:theme xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" name="PPIL4">'
            '<a:themeElements><a:clrScheme name="PPIL4">' + scheme + '</a:clrScheme>'
            '<a:fontScheme name="PPIL4">' + font + '</a:fontScheme>' + fmt
            + '</a:themeElements><a:objectDefaults/><a:extraClrSchemeLst/></a:theme>')


def export_pptx():
    out = os.path.join(HERE, "PPIL4_deck.pptx")
    imgs = render_slides()
    n = len(imgs)

    ct = [XML, '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
          '<Default Extension="xml" ContentType="application/xml"/>',
          '<Default Extension="png" ContentType="image/png"/>']
    for part, kind in [("/ppt/presentation.xml", "presentationml.presentation.main"),
                       ("/ppt/slideMasters/slideMaster1.xml", "presentationml.slideMaster"),
                       ("/ppt/slideLayouts/slideLayout1.xml", "presentationml.slideLayout"),
                       ("/docProps/app.xml", "extended-properties")]:
        ct.append('<Override PartName="%s" ContentType="application/vnd.openxmlformats-officedocument.%s+xml"/>' % (part, kind))
    ct.append('<Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/>')
    ct.append('<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>')
    ct += ['<Override PartName="/ppt/slides/slide%d.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>' % (i + 1)
           for i in range(n)]
    ct.append('</Types>')

    pres = (XML + '<p:presentation %s saveSubsetFonts="1">' % NS
            + '<p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId1"/></p:sldMasterIdLst><p:sldIdLst>'
            + "".join('<p:sldId id="%d" r:id="rId%d"/>' % (256 + i, i + 2) for i in range(n))
            + '</p:sldIdLst><p:sldSz cx="%d" cy="%d"/><p:notesSz cx="%d" cy="%d"/></p:presentation>' % (W, H, H, W))

    pres_rels = [("rId1", REL + "slideMaster", "slideMasters/slideMaster1.xml")]
    pres_rels += [("rId%d" % (i + 2), REL + "slide", "slides/slide%d.xml" % (i + 1)) for i in range(n)]
    pres_rels += [("rId%d" % (n + 2), REL + "theme", "theme/theme1.xml")]

    master = (XML + '<p:sldMaster %s><p:cSld>%s</p:cSld>' % (NS, EMPTY_TREE)
              + '<p:clrMap bg1="lt1" tx1="dk1" bg2="lt2" tx2="dk2" accent1="accent1" accent2="accent2" '
                'accent3="accent3" accent4="accent4" accent5="accent5" accent6="accent6" hlink="hlink" '
                'folHlink="folHlink"/><p:sldLayoutIdLst><p:sldLayoutId id="2147483649" r:id="rId1"/>'
                '</p:sldLayoutIdLst></p:sldMaster>')
    layout = (XML + '<p:sldLayout %s type="blank" preserve="1"><p:cSld name="Blank">%s</p:cSld>' % (NS, EMPTY_TREE)
              + '<p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sldLayout>')

    def slide_xml(i):
        return (XML + '<p:sld %s><p:cSld><p:spTree>' % NS
                + '<p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>'
                  '<p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/>'
                  '<a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>'
                + '<p:pic><p:nvPicPr><p:cNvPr id="2" name="Slide %d"/>' % (i + 1)
                + '<p:cNvPicPr><a:picLocks noChangeAspect="1"/></p:cNvPicPr><p:nvPr/></p:nvPicPr>'
                  '<p:blipFill><a:blip r:embed="rId1"/><a:stretch><a:fillRect/></a:stretch></p:blipFill>'
                + '<p:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="%d" cy="%d"/></a:xfrm>' % (W, H)
                + '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic>'
                  '</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>')

    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    core = (XML + '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
            'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            '<dc:title>能不能通过建一座桥治疗胰腺癌?</dc:title>'
            '<dcterms:created xsi:type="dcterms:W3CDTF">%s</dcterms:created>'
            '<dcterms:modified xsi:type="dcterms:W3CDTF">%s</dcterms:modified></cp:coreProperties>' % (now, now))
    app = (XML + '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
           'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
           '<Slides>%d</Slides><Application>deck/export.py</Application></Properties>' % n)

    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", "".join(ct))
        z.writestr("_rels/.rels", _rels([
            ("rId1", REL + "officeDocument", "ppt/presentation.xml"),
            ("rId2", "http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties", "docProps/core.xml"),
            ("rId3", REL + "extended-properties", "docProps/app.xml")]))
        z.writestr("docProps/core.xml", core)
        z.writestr("docProps/app.xml", app)
        z.writestr("ppt/presentation.xml", pres)
        z.writestr("ppt/_rels/presentation.xml.rels", _rels(pres_rels))
        z.writestr("ppt/theme/theme1.xml", _theme())
        z.writestr("ppt/slideMasters/slideMaster1.xml", master)
        z.writestr("ppt/slideMasters/_rels/slideMaster1.xml.rels", _rels([
            ("rId1", REL + "slideLayout", "../slideLayouts/slideLayout1.xml"),
            ("rId2", REL + "theme", "../theme/theme1.xml")]))
        z.writestr("ppt/slideLayouts/slideLayout1.xml", layout)
        z.writestr("ppt/slideLayouts/_rels/slideLayout1.xml.rels", _rels([
            ("rId1", REL + "slideMaster", "../slideMasters/slideMaster1.xml")]))
        for i, img in enumerate(imgs):
            z.write(img, "ppt/media/image%d.png" % (i + 1))
            z.writestr("ppt/slides/slide%d.xml" % (i + 1), slide_xml(i))
            z.writestr("ppt/slides/_rels/slide%d.xml.rels" % (i + 1), _rels([
                ("rId1", REL + "image", "../media/image%d.png" % (i + 1))]))
    print("  PPTX  %s  (%.1f MB, %d slides)" % (out, os.path.getsize(out) / 1048576, n))
    return out


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "both"
    try:
        if what in ("pdf", "both"):
            export_pdf()
        if what in ("pptx", "both"):
            export_pptx()
    finally:
        shutil.rmtree(TMP, ignore_errors=True)
