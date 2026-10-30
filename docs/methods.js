// The PCA panels: a hover label like the plots on the Genes page, the hovered sample enlarged,
// its patient's other sample and their line picked out, and a legend hover that dims the other groups.
const tip = document.createElement("div");
tip.className = "tip"; tip.hidden = true; document.body.appendChild(tip);

function place(e) {
  tip.style.left = `${e.pageX + 14}px`;
  tip.style.top = `${e.pageY - 12}px`;
}
document.querySelectorAll(".pca-panel").forEach((panel) => {
  const circles = [...panel.querySelectorAll("circle")], lines = [...panel.querySelectorAll("line.pair")];
  const show = (c, e) => {
    tip.innerHTML = `<b>${c.dataset.label}</b><br><span class="swatch" style="background:${c.dataset.color}"></span>${c.dataset.group}<br><span class="dim">${c.dataset.pc}</span>`;
    tip.hidden = false;
    if (e) place(e); else { const r = c.getBoundingClientRect(); tip.style.left = `${r.right + window.scrollX + 8}px`; tip.style.top = `${r.top + window.scrollY - 12}px`; }
    circles.forEach((o) => o.classList.toggle("mate", o !== c && !!c.dataset.patient && o.dataset.patient === c.dataset.patient));
    lines.forEach((l) => l.classList.toggle("lit", !!c.dataset.patient && l.dataset.patient === c.dataset.patient));
    c.classList.add("hot");
  };
  const hide = () => {
    tip.hidden = true;
    circles.forEach((o) => o.classList.remove("hot", "mate"));
    lines.forEach((l) => l.classList.remove("lit"));
  };
  circles.forEach((c) => {
    c.addEventListener("mouseenter", (e) => show(c, e));
    c.addEventListener("mousemove", place);
    c.addEventListener("mouseleave", hide);
    c.addEventListener("focus", () => show(c));
    c.addEventListener("blur", hide);
  });
  panel.querySelectorAll(".legend [data-group]").forEach((el) => {
    el.addEventListener("mouseenter", () => circles.forEach((c) => c.classList.toggle("faded", c.dataset.group !== el.dataset.group)));
    el.addEventListener("mouseleave", () => circles.forEach((c) => c.classList.remove("faded")));
  });
});
