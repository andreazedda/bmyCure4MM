"use strict";
const confirmation = document.querySelector("[data-token-form]");
if (confirmation) {
  const token = new URLSearchParams(location.hash.slice(1)).get("token");
  history.replaceState(null, "", location.pathname);
  if (token) confirmation.querySelector('[name="token"]').value = token;
}
