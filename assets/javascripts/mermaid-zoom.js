// Click-to-zoom for Mermaid SVG diagrams.
// Opens a fullscreen overlay with wheel zoom + drag pan.

function setupMermaidZoom() {
  const svgs = Array.from(document.querySelectorAll(".mermaid svg")).filter(
    (svg) => !svg.getAttribute("data-zoom-bound"),
  );
  if (!svgs.length) return;

  function openOverlay(svg) {
    const overlay = document.createElement("div");
    overlay.className = "mermaid-zoom-overlay";
    overlay.innerHTML = `
      <div class="mermaid-zoom-toolbar">
        <button type="button" class="mermaid-zoom-btn" data-action="reset">Reset</button>
        <button type="button" class="mermaid-zoom-btn" data-action="close">Close</button>
      </div>
      <div class="mermaid-zoom-stage">
        <div class="mermaid-zoom-canvas" tabindex="0" role="dialog" aria-label="Mermaid diagram zoom">
          ${svg.outerHTML}
        </div>
      </div>
    `;

    document.body.appendChild(overlay);

    const canvas = overlay.querySelector(".mermaid-zoom-canvas");
    const zoomSvg = canvas.querySelector("svg");

    // Ensure the cloned SVG doesn't inherit our binding flag.
    zoomSvg.removeAttribute("data-zoom-bound");

    let scale = 1;
    let tx = 0;
    let ty = 0;
    let dragging = false;
    let lastX = 0;
    let lastY = 0;

    function apply() {
      canvas.style.transform = `translate(${tx}px, ${ty}px) scale(${scale})`;
    }

    function reset() {
      scale = 1;
      tx = 0;
      ty = 0;
      apply();
    }

    function close() {
      overlay.remove();
    }

    overlay.addEventListener("click", (e) => {
      if (e.target === overlay) close();
    });

    overlay.querySelector("[data-action='close']").addEventListener("click", close);
    overlay.querySelector("[data-action='reset']").addEventListener("click", reset);

    window.addEventListener(
      "keydown",
      (e) => {
        if (e.key === "Escape") close();
      },
      { once: true },
    );

    canvas.addEventListener("wheel", (e) => {
      e.preventDefault();
      const delta = Math.sign(e.deltaY);
      const factor = delta > 0 ? 0.9 : 1.1;
      scale = Math.max(0.2, Math.min(6, scale * factor));
      apply();
    });

    canvas.addEventListener("pointerdown", (e) => {
      dragging = true;
      lastX = e.clientX;
      lastY = e.clientY;
      canvas.setPointerCapture(e.pointerId);
    });

    canvas.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      tx += e.clientX - lastX;
      ty += e.clientY - lastY;
      lastX = e.clientX;
      lastY = e.clientY;
      apply();
    });

    canvas.addEventListener("pointerup", () => {
      dragging = false;
    });

    canvas.addEventListener("dblclick", reset);

    // Prevent the zoomed SVG from being constrained.
    zoomSvg.style.maxWidth = "none";
    zoomSvg.style.maxHeight = "none";

    reset();
  }

  for (const svg of svgs) {
    svg.setAttribute("data-zoom-bound", "true");
    svg.style.cursor = "zoom-in";
    svg.addEventListener("click", (e) => {
      // If user selects text or clicks links inside SVG, don't hijack.
      if (e.defaultPrevented) return;
      e.preventDefault();
      openOverlay(svg);
    });
  }
}

if (typeof document$ !== "undefined" && document$?.subscribe) {
  document$.subscribe(() => setupMermaidZoom());
}

setupMermaidZoom();

