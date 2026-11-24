document.addEventListener("DOMContentLoaded", () => {
  // Theme switcher
  const picker = document.getElementById("theme-picker");
  const saved = localStorage.getItem("uw-theme") || "light";
  document.documentElement.setAttribute("data-theme", saved);
  if (picker) {
    picker.value = saved;
    picker.addEventListener("change", (e) => {
      const val = e.target.value;
      document.documentElement.setAttribute("data-theme", val);
      localStorage.setItem("uw-theme", val);
    });
  }

  // Copy buttons on code blocks
  document.querySelectorAll("pre").forEach((pre) => {
    const btn = document.createElement("button");
    btn.className = "copy-btn";
    btn.textContent = "Copy";
    btn.addEventListener("click", () => {
      const code = pre.querySelector("code") ? pre.querySelector("code").innerText : pre.innerText;
      navigator.clipboard.writeText(code).then(() => {
        btn.textContent = "Copied!";
        setTimeout(() => { btn.textContent = "Copy"; }, 1400);
      });
    });
    pre.appendChild(btn);
  });

  // Sidebar filter search
  const search = document.getElementById("sidebar-search");
  if (search) {
    search.addEventListener("input", (e) => {
      const q = e.target.value.toLowerCase().trim();
      document.querySelectorAll(".nav-list li").forEach((li) => {
        const txt = li.textContent.toLowerCase();
        li.style.display = (!q || txt.includes(q)) ? "" : "none";
      });
    });
  }

  // Render KaTeX math if loaded
  if (window.renderMathInElement) {
    window.renderMathInElement(document.body, {
      delimiters: [
        { left: "$$", right: "$$", display: true },
        { left: "\\[", right: "\\]", display: true },
        { left: "\\(", right: "\\)", display: false },
        { left: "$", right: "$", display: false }
      ],
      throwOnError: false
    });
  }
});
