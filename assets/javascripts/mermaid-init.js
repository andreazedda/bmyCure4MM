/* Mermaid bootstrap for MkDocs Material (supports instant navigation). */

function renderMermaid() {
  if (!window.mermaid) return;

  const scheme = document.body.getAttribute("data-md-color-scheme");
  const isDark = scheme === "slate";

  window.mermaid.initialize({
    startOnLoad: false,
    theme: isDark ? "dark" : "default",
    flowchart: {
      nodeSpacing: 30,
      rankSpacing: 30,
    },
    themeVariables: {
      fontSize: "16px",
    },
  });

  const nodes = Array.from(document.querySelectorAll(".mermaid")).filter((node) => {
    if (node.getAttribute("data-processed")) return false;
    // If mermaid already replaced it with SVG, don't touch it again.
    if (node.querySelector?.("svg")) return false;
    return true;
  });
  if (!nodes.length) return;

  nodes.forEach((node, idx) => {
    const id = `mermaid-${Date.now()}-${idx}`;
    const source = (node.textContent || "").trim();
    if (!source) return;

    window.mermaid
      .render(id, source)
      .then(({ svg }) => {
        node.innerHTML = svg;
        node.setAttribute("data-processed", "true");
      })
      .catch((err) => {
        node.setAttribute("data-processed", "true");
        // eslint-disable-next-line no-console
        console.error("Mermaid render failed:", err);
      });
  });
}

if (typeof document$ !== "undefined" && document$?.subscribe) {
  document$.subscribe(() => renderMermaid());
}

renderMermaid();
