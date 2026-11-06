"use strict";
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
    } catch {
      button.textContent = "Select code / Seleziona il codice";
    }
  });
  container.append(button);
});
