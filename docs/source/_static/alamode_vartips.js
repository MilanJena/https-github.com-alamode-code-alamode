// Hover/focus pop-ups for links to input-variable entries, and hover pop-ups for
// tags and &field names in ALAMODE input code blocks (alamode_vartips.py).
// The entries live in _static/alamode_vartips_data.js, loaded on first use.
(function () {
  "use strict";
  const ENTRY = /#(alm|anphon)-[\w-]+$/;
  const CODE = "span.na-Alamode, span.kn-Alamode"; // tags, &fields (alamode_lexer.py)
  const ALM_MODES = ["suggest", "optimize", "opt"];
  const root = new URL(document.documentElement.dataset.content_root || "./", location.href).href;
  let tip = null, anchor = null, showTimer = 0, hideTimer = 0, loading = null, mathjax = null;
  const blocks = new WeakMap(); // <pre> -> {prog, qha}

  function loadScript(src) {
    return new Promise(function (resolve) {
      const s = document.createElement("script");
      s.src = src;
      s.onload = s.onerror = resolve;
      document.head.appendChild(s);
    });
  }

  function load() {
    if (window.ALAMODE_VARTIPS) return Promise.resolve();
    if (!loading) loading = loadScript(root + "_static/alamode_vartips_data.js");
    return loading;
  }

  // "anphondir/inputanphon.html#anphon-ismear" for a link to an entry on another page
  function linkKey(a) {
    if (!a.classList.contains("internal")) return null;
    const u = new URL(a.getAttribute("href"), location.href);
    if (!ENTRY.test(u.hash) || u.pathname === location.pathname || !u.href.startsWith(root)) return null;
    return u.href.slice(root.length);
  }

  // Program of an input block: fields only one program has, then MODE, then tags
  // only one program has; null if still undecided.
  function classify(pre) {
    const V = window.ALAMODE_VARTAGS;
    const only = function (kind, names) {
      const a = names.some(function (n) { return n in V.alm[kind] && !(n in V.anphon[kind]); });
      const b = names.some(function (n) { return n in V.anphon[kind] && !(n in V.alm[kind]); });
      return a === b ? null : a ? "alm" : "anphon";
    };
    const text = function (sel) {
      return Array.from(pre.querySelectorAll(sel), function (s) { return s.textContent; });
    };
    const fields = text("span.kn-Alamode").map(function (f) { return f.slice(1).toLowerCase(); });
    const mode = /(?:^|[\s;])MODE\s*=\s*([^\s;#]+)/m.exec(pre.textContent.replace(/#.*$/gm, ""));
    const prog = only("fields", fields) ||
      (mode ? (ALM_MODES.includes(mode[1].toLowerCase()) ? "alm" : "anphon") : null) ||
      only("tags", text("span.na-Alamode"));
    return { prog: prog, qha: fields.includes("qha") };
  }

  function codeKey(span) {
    const V = window.ALAMODE_VARTAGS, pre = span.closest("pre");
    if (!V || !pre) return null;
    if (!blocks.has(pre)) blocks.set(pre, classify(pre));
    const b = blocks.get(pre), field = span.classList.contains("kn-Alamode");
    const name = field ? span.textContent.slice(1).toLowerCase() : span.textContent;
    const lookup = function (prog) {
      const m = V[prog];
      if (field) return m.fields[name];
      return (prog === "anphon" && b.qha && m.alt[name]) || m.tags[name];
    };
    if (b.prog) return lookup(b.prog) || null;
    // undecided block: only names that one program alone knows
    const kind = field ? "fields" : "tags";
    const inAlm = name in V.alm[kind], inAnphon = name in V.anphon[kind];
    return inAlm === inAnphon ? null : lookup(inAlm ? "alm" : "anphon") || null;
  }

  function hide() {
    clearTimeout(showTimer);
    if (tip) tip.hidden = true;
    anchor = null;
  }

  // Typeset the tip's math, loading MathJax first on pages that have none.
  function typeset() {
    if (!tip.querySelector(".math")) return;
    if (!mathjax && !(window.MathJax && MathJax.startup)) {
      const cfg = window.ALAMODE_MATHJAX || {};
      if (!cfg.path) return;
      const src = /^(\w+:)?\//.test(cfg.path) ? cfg.path : root + "_static/" + cfg.path;
      // The page's own MathJax (deferred, after this script) may still be loading.
      const own = document.querySelector('script[src="' + cfg.path + '"]');
      if (own) {
        mathjax = new Promise(function (resolve) { own.addEventListener("load", resolve); });
      } else {
        window.MathJax = Object.assign({}, cfg.config || {}, { startup: { typeset: false } });
        mathjax = loadScript(src);
      }
      mathjax = mathjax.then(function () {
        return window.MathJax.startup && MathJax.startup.promise;
      });
    }
    const target = tip;
    Promise.resolve(mathjax).then(function () {
      if (window.MathJax && MathJax.typesetPromise) MathJax.typesetPromise([target]);
    });
  }

  function show(el, key) {
    const html = (window.ALAMODE_VARTIPS || {})[key];
    if (!html || anchor !== el) return;
    if (!tip) {
      tip = document.createElement("div");
      tip.className = "alamode-vartip";
      tip.setAttribute("role", "tooltip");
      tip.addEventListener("mouseenter", function () { clearTimeout(hideTimer); });
      tip.addEventListener("mouseleave", scheduleHide);
      document.body.appendChild(tip);
    }
    if (window.MathJax && MathJax.typesetClear) MathJax.typesetClear([tip]);
    tip.innerHTML = html;
    tip.querySelectorAll("[id]").forEach(function (e) { e.removeAttribute("id"); });
    // hrefs are relative to the input page: make them absolute
    const base = new URL(key.split("#")[0], root);
    tip.querySelectorAll("[href]").forEach(function (e) {
      e.setAttribute("href", new URL(e.getAttribute("href"), base).href);
    });
    tip.querySelectorAll("[src]").forEach(function (e) {
      e.setAttribute("src", new URL(e.getAttribute("src"), base).href);
    });
    tip.hidden = false;
    tip.scrollTop = 0;
    const r = el.getBoundingClientRect(), w = tip.offsetWidth, h = tip.offsetHeight;
    const left = Math.max(8, Math.min(r.left, document.documentElement.clientWidth - w - 8));
    const top = r.bottom + 6 + h < window.innerHeight ? r.bottom + 6 : Math.max(8, r.top - 6 - h);
    tip.style.left = left + window.scrollX + "px";
    tip.style.top = top + window.scrollY + "px";
    typeset();
  }

  function scheduleHide() {
    clearTimeout(showTimer);
    clearTimeout(hideTimer);
    hideTimer = setTimeout(hide, 300);
  }

  function target(e) {
    return e.target.closest ? e.target.closest("a.reference, " + CODE) : null;
  }

  function enter(e) {
    const el = target(e);
    if (!el) return;
    if (el === anchor) {
      clearTimeout(hideTimer);
      if (tip && !tip.hidden) return; // already showing
    }
    const isLink = el.tagName === "A";
    if (isLink ? !linkKey(el) : e.type === "focusin") return;
    clearTimeout(hideTimer);
    clearTimeout(showTimer);
    anchor = el;
    showTimer = setTimeout(function () {
      load().then(function () {
        const key = isLink ? linkKey(el) : codeKey(el);
        if (!key) {
          if (anchor === el) hide();
          return;
        }
        if (!isLink) el.classList.add("alamode-vartip-target");
        show(el, key);
      });
    }, e.type === "focusin" ? 0 : 250);
  }

  function leave(e) {
    const el = target(e);
    if (el && el === anchor && !(e.relatedTarget && el.contains(e.relatedTarget))) scheduleHide();
  }

  document.addEventListener("mouseover", enter);
  document.addEventListener("focusin", enter);
  document.addEventListener("mouseout", leave);
  document.addEventListener("focusout", leave);
  document.addEventListener("keydown", function (e) {
    if (e.key === "Escape") hide();
  });
})();
