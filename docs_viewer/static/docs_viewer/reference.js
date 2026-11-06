"use strict";
// Keep reading content ahead of the long contents list on narrow screens.
document.querySelectorAll('.doc-navigation').forEach((navigation) => {
  const breakpoint = window.matchMedia('(min-width: 768px)');
  navigation.open = breakpoint.matches;
  breakpoint.addEventListener('change', (event) => { navigation.open = event.matches; });
});
document.querySelectorAll(".doc-content .highlight").forEach((container) => {
  const code = container.querySelector("pre");
  if (!code) return;
  const button = document.createElement("button");
  button.type = "button";
  button.className = "doc-copy";
  button.innerHTML = '<span class="t-en">Copy code</span><span class="t-it">Copia codice</span>';
  button.addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText(code.textContent);
      const italian = document.documentElement.lang.startsWith("it");
      button.textContent = italian ? "Copiato" : "Copied";
      if (window.updateLive) window.updateLive(button.textContent);
    } catch {
      const selection = window.getSelection();
      const range = document.createRange();
      range.selectNodeContents(code);
      selection.removeAllRanges();
      selection.addRange(range);
      button.textContent = document.documentElement.lang.startsWith('it') ? 'Codice selezionato: usa Copia' : 'Code selected: use Copy';
      if (window.updateLive) window.updateLive(button.textContent);
    }
  });
  container.append(button);
});

// Keep table semantics while giving the scroll region a keyboard target.
document.querySelectorAll('.doc-content table').forEach((table) => {
  const scroller = document.createElement('div');
  scroller.className = 'doc-table-scroll';
  table.before(scroller);
  scroller.append(table);
});
// Native arrow-key scrolling remains available for overflowing references.
const updateDocScrollRegions = () => {
  document.querySelectorAll('.doc-content .doc-table-scroll, .doc-content pre, .doc-content .arithmatex, .doc-content .mermaid').forEach((node) => {
    if (node.scrollWidth > node.clientWidth + 1) {
      node.tabIndex = 0;
      node.setAttribute('role', 'region');
      node.setAttribute('aria-label', document.documentElement.lang.startsWith('it') ? 'Contenuto scorrevole orizzontalmente' : 'Horizontally scrollable content');
    } else {
      node.removeAttribute('tabindex');
      node.removeAttribute('role');
      node.removeAttribute('aria-label');
    }
  });
};
window.addEventListener('resize', updateDocScrollRegions);
window.addEventListener('load', updateDocScrollRegions);
document.addEventListener('language-changed', updateDocScrollRegions);
updateDocScrollRegions();

const article = document.querySelector('.doc-content');
if (article) new MutationObserver(updateDocScrollRegions).observe(article, {childList: true, subtree: true});
