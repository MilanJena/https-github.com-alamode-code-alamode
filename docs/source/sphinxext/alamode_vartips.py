"""Links to input variables and hover pop-ups with their reference entries.

``:alm:`NORDER``` and ``:anphon:`ISMEAR``` render as inline literals linking to
the tag's entry in almdir/inputalm.rst or anphondir/inputanphon.rst
(``:anphon:`text <TAG or label>``` sets the text). The tag -> label map is read
from the bullet line under each ``.. _alm_*:`` / ``.. _anphon_*:`` label, so
entries covering several tags (``TMIN, TMAX, DT-tags``) work; an unknown tag
warns.

HTML builders also write ``_static/alamode_vartips_data.js`` with the rendered
entries, keyed by ``page#anchor`` (``*_field`` section labels get the title and
first paragraph), plus per-program tag/field maps and the MathJax URL;
``_static/alamode_vartips.js`` shows the entries when hovering any link to an
entry, or a tag or ``&field`` name in an ALAMODE input code block. ``python alamode_vartips.py`` runs a
self-check on the two input pages.
"""

import json
import re
from pathlib import Path

from docutils import nodes
from docutils.utils import new_document
from sphinx import addnodes
from sphinx.util import logging
from sphinx.util.nodes import make_refnode, split_explicit_title

logger = logging.getLogger(__name__)

SRCDIR = Path(__file__).resolve().parent.parent
PAGES = {"alm": "almdir/inputalm", "anphon": "anphondir/inputanphon"}
LABEL_RE = re.compile(r"^\.\. _((?:alm|anphon)_\w+):\s*$")
TAG_RE = re.compile(r"^[A-Z][A-Z0-9_]+$")
FIELD_RE = re.compile(r'^"&(\w+)"-field')


def parse_tags(path):
    """rst file -> ({TAG: first label of its entry}, set of all labels,
    {field: label on its section or None}, {TAG: label of a second entry})."""
    tags, labels, fields, alt, pending = {}, set(), {}, {}, []
    for line in Path(path).read_text().splitlines():
        m = LABEL_RE.match(line)
        if m:
            pending.append(m.group(1))
            labels.add(m.group(1))
            continue
        f = FIELD_RE.match(line)
        if f:
            fields[f.group(1)] = pending[0] if pending else None
            pending = []
        elif line.strip() and pending:
            head = re.split(r"-tag| =| :", line.replace("**", "").lstrip("* "))[0]
            for t in head.split(","):
                t = t.strip()
                if TAG_RE.match(t):
                    if t in tags:
                        alt.setdefault(t, pending[0])
                    tags.setdefault(t, pending[0])
            pending = []
    return tags, labels, fields, alt


_maps = {}


def tag_map(prog):
    if prog not in _maps:
        _maps[prog] = parse_tags(SRCDIR / f"{PAGES[prog]}.rst")
    return _maps[prog]


def make_role(prog):
    def role(name, rawtext, text, lineno, inliner, options=None, content=None):
        env = inliner.document.settings.env
        _, title, target = split_explicit_title(text)
        # rebuild this page when the input reference changes (labels may move)
        env.note_dependency(str(SRCDIR / f"{PAGES[prog]}.rst"))
        tags, labels = tag_map(prog)[:2]
        label = target if target in labels else tags.get(target.upper())
        if label is None:
            logger.warning(
                "unknown %s input variable: %r",
                prog,
                target,
                location=(env.docname, lineno),
            )
            return [nodes.literal(rawtext, title)], []
        node = addnodes.pending_xref(
            rawtext,
            nodes.literal(rawtext, title, classes=["xref", "std", "std-ref"]),
            refdomain="",
            reftype="vartip",
            reftarget=label,
            refexplicit=True,
            refwarn=True,
            refdoc=env.docname,
        )
        return [node], []

    return role


def resolve(app, env, node, contnode):
    """Link to the label like :ref: does, but keep the literal (:ref: replaces
    an explicit title by an inline; std "keyword" only knows section labels)."""
    if node["reftype"] != "vartip":
        return None
    try:
        docname, labelid = env.domains.standard_domain.anonlabels[node["reftarget"]]
    except KeyError:
        return None  # refwarn: Sphinx reports the undefined label
    return make_refnode(app.builder, node["refdoc"], docname, labelid, contnode)


def _entry_ids(n):
    ids = n.get("ids", ()) if isinstance(n, nodes.Element) else ()
    return [i for i in ids if i.startswith(("alm-", "anphon-"))]


def _render(builder, group):
    """Resolved doctree nodes -> HTML (render_partial's transforms choke on
    Sphinx admonitions, and the nodes need none)."""
    doc = new_document("<vartip>", builder._settings)
    for m in group:
        doc.append(m.deepcopy())
    visitor = builder.create_translator(doc, builder)
    doc.walkabout(visitor)
    return "".join(visitor.fragment)


def _section_tip(sec):
    """Section -> [its title in bold, its first paragraph]."""
    head = nodes.paragraph("", "", nodes.strong("", "", *sec[0].deepcopy().children))
    para = next(iter(sec.findall(nodes.paragraph)), None)
    return [head] if para is None else [head, para]


def write_tips(app, exception):
    builder = app.builder
    # singlehtml targets are "#document-..." fragments; the JS keys assume pages
    if (
        exception is not None
        or builder.format != "html"
        or builder.name == "singlehtml"
    ):
        return
    tips = {}
    for docname in PAGES.values():
        if docname not in app.env.found_docs:
            continue
        doctree = app.env.get_and_resolve_doctree(docname, builder)
        builder.current_docname = docname
        page = builder.get_target_uri(docname)
        for parent in doctree.findall(nodes.Element):
            kids = parent.children
            for i, n in enumerate(kids):
                ids = _entry_ids(n)
                if isinstance(n, nodes.section):  # only *_field labels
                    ids = [x for x in ids if x.endswith("-field")]
                if not ids:
                    continue
                if isinstance(n, nodes.section):
                    group = _section_tip(n)
                else:
                    group = [n]
                    for m in kids[i + 1 :]:
                        if isinstance(
                            m, (nodes.transition, nodes.section)
                        ) or _entry_ids(m):
                            break
                        group.append(m)
                html = _render(builder, group)
                for id_ in ids:
                    tips[f"{page}#{id_}"] = html

    # Tag and field maps for the pop-ups on ALAMODE input code blocks.
    anon = app.env.domains.standard_domain.anonlabels

    def key(label):
        if label not in anon:
            return None
        doc, id_ = anon[label]
        return f"{builder.get_target_uri(doc)}#{id_}"

    vartags = {}
    for prog in PAGES:
        tags, _, fields, alt = tag_map(prog)
        vartags[prog] = {
            "tags": {t: key(lb) for t, lb in tags.items() if key(lb)},
            "alt": {t: key(lb) for t, lb in alt.items() if key(lb)},
            "fields": {f: key(lb) if lb else "" for f, lb in fields.items()},
        }
    mathjax = {"path": app.config.mathjax_path}
    mathjax["config"] = app.config.mathjax4_config or app.config.mathjax3_config
    out = Path(builder.outdir) / "_static" / "alamode_vartips_data.js"
    out.write_text(
        "window.ALAMODE_VARTIPS = "
        + json.dumps(tips, ensure_ascii=False)
        + ";\nwindow.ALAMODE_VARTAGS = "
        + json.dumps(vartags, ensure_ascii=False)
        + ";\nwindow.ALAMODE_MATHJAX = "
        + json.dumps(mathjax)
        + ";\n",
        encoding="utf-8",
    )


def setup(app):
    for prog in PAGES:
        app.add_role(prog, make_role(prog))
    app.add_js_file("alamode_vartips.js", loading_method="defer")
    app.add_css_file("alamode_vartips.css")
    app.connect("builder-inited", lambda app: _maps.clear())  # re-parse per build
    app.connect("missing-reference", resolve)
    app.connect("build-finished", write_tips)
    return {"parallel_read_safe": True, "parallel_write_safe": True}


if __name__ == "__main__":
    for prog in PAGES:
        tags, labels, fields, alt = tag_map(prog)
        unmapped = labels - set(tags.values()) - set(alt.values())
        unmapped -= set(fields.values())
        print(
            prog, len(tags), "tags,", len(labels), "labels; no tag:", sorted(unmapped)
        )
        print("  fields:", fields, "\n  second entries:", alt)
    a, _, f, alt = tag_map("anphon")
    assert a["TMAX"] == a["DT"] == "anphon_tmin" and a["DELTA_E"] == "anphon_emin"
    assert a["ANIME_FORMAT"] == "anphon_anime_format" and a["ISMEAR"] == "anphon_ismear"
    assert a["RENORM_3TO2ND"] == "anphon_elastic_const"
    assert a["RELAX_STR"] == "anphon_relax_str"
    assert alt["RELAX_STR"] == "anphon_qha_relax_str"
    assert f["cell"] == "anphon_cell_field" and f["kappa"] is None
    b, _, f, _ = tag_map("alm")
    assert (
        b["PERIODIC_IMAGE_CONV"] == "alm_periodic_image_conv"
        and b["NEND"] == "alm_nstart"
    )
    assert "interaction" in f and "kpoint" not in f
    print("alamode_vartips self-check OK")
