// MathJax config for MkDocs Material + instant navigation.
// Uses pymdownx.arithmatex (generic mode) to wrap math spans.

window.MathJax = {
  tex: {
    inlineMath: [["$", "$"], ["\\(", "\\)"]],
    displayMath: [["$$", "$$"], ["\\[", "\\]"]],
    processEscapes: true,
  },
  options: {
    skipHtmlTags: ["script", "noscript", "style", "textarea", "pre", "code"],
    ignoreHtmlClass: ".*",
    processHtmlClass: "arithmatex",
  },
};

let pendingTypeset = true;

function typesetMath() {
  if (!window.MathJax?.typesetPromise) {
    pendingTypeset = true;
    return;
  }
  pendingTypeset = false;
  window.MathJax.typesetPromise();
}

if (typeof document$ !== "undefined" && document$?.subscribe) {
  document$.subscribe(() => typesetMath());
}

// First load (or if user navigates before MathJax finishes loading).
window.addEventListener("load", () => {
  if (pendingTypeset) typesetMath();
});
