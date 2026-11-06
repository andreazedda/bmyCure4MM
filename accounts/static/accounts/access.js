"use strict";
const confirmation = document.querySelector("[data-token-form]");
if (confirmation) {
  const token = new URLSearchParams(location.hash.slice(1)).get("token");
  history.replaceState(null, "", location.pathname);
  if (token) confirmation.querySelector('[name="token"]').value = token;
}

const localizeAccountChoices = () => {
  const lang = document.documentElement.lang === "it" ? "it" : "en";
  document.querySelectorAll("option[data-option-en]").forEach(option => {
    option.textContent = option.getAttribute("data-option-" + lang);
  });
};
localizeAccountChoices();
document.addEventListener("language-changed", localizeAccountChoices);
const errorSummary = document.querySelector("[data-error-summary]");
if (errorSummary) errorSummary.focus();
if (confirmation && !confirmation.querySelector('[name="token"]').value) {
  const missing = document.querySelector("[data-token-missing]");
  if (missing) missing.hidden = false;
  confirmation.querySelector('button[type="submit"]').disabled = true;
}
