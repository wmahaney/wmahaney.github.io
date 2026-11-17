// The side panel on the home page that holds the isogeny graphs (index.qmd).
// The tab on the right edge opens it; the × button, the Escape key, or the tab again
// close it. On a wide screen the page moves over to make room; on a narrow one the
// panel covers it. While closed the panel is "inert" (nothing in it can be focused),
// and the graph pauses because it is off screen (see isogeny-hero.js).

const panel = document.getElementById("wm-side");
const tab = document.querySelector(".wm-side-tab");

if (panel && tab) {
  const close = panel.querySelector(".wm-side-close");
  // Quarto puts the page inside containers that start their own stacking order, which
  // would leave the menu bar on top of the panel; at the end of <body> it is on top.
  document.body.append(tab, panel);

  const setOpen = (open, focus = true) => {
    document.body.classList.toggle("wm-side-open", open);
    tab.setAttribute("aria-expanded", String(open));
    panel.inert = !open;
    if (focus) (open ? close : tab).focus({ preventScroll: true });
  };

  setOpen(false, false);
  tab.hidden = false;
  document.body.classList.add("wm-side-ready");

  tab.addEventListener("click", () => setOpen(!document.body.classList.contains("wm-side-open")));
  close.addEventListener("click", () => setOpen(false));
  // Where the panel covers the page, a click on the page closes it.
  const covers = window.matchMedia("(max-width: 1599.98px)");
  document.addEventListener("click", (e) => {
    if (covers.matches && document.body.classList.contains("wm-side-open")
        && !panel.contains(e.target) && !tab.contains(e.target)) setOpen(false, false);
  });
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && document.body.classList.contains("wm-side-open")
        && !panel.querySelector(".wm-graph-pop:not([hidden])")) setOpen(false);
  });
}
