(() => {
  const header = document.querySelector('.site-tagline');
  const data = document.getElementById('header-lines');
  if (!header || !data) return;
  try {
    const lines = JSON.parse(data.textContent);
    if (!Array.isArray(lines) || !lines.length) return;
    const choice = lines[Math.floor(Math.random() * lines.length)];
    if (!choice || typeof choice.text !== 'string') return;
    header.textContent = choice.text;
    header.dataset.headerId = choice.id;
  } catch (_) { /* Keep the static header if data is unavailable. */ }
})();
