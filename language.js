(() => {
  const key = 'du-pei-chinese-script';
  function preferredScript() {
    try {
      const saved = localStorage.getItem(key);
      if (saved === 'Hans' || saved === 'Hant') return saved;
    } catch (_) { /* Storage may be disabled; browser language still works. */ }
    for (const language of navigator.languages || [navigator.language]) {
      const tag = (language || '').toLowerCase();
      if (!/^zh(?:-|$)/.test(tag)) continue;
      if (tag.includes('-hant')) return 'Hant';
      if (tag.includes('-hans')) return 'Hans';
      return /-(tw|hk|mo)(?:-|$)/.test(tag) ? 'Hant' : 'Hans';
    }
    return 'Hans';
  }
  if (document.documentElement.lang === 'en') {
    const chineseTab = document.querySelector('nav a[lang="zh"]');
    if (chineseTab) chineseTab.href = preferredScript() === 'Hant' ? 'zh-hant.html' : 'zh.html';
  }
  document.querySelectorAll('.language-switch a').forEach(link => {
    link.addEventListener('click', () => {
      try { localStorage.setItem(key, link.lang === 'zh-Hant' ? 'Hant' : 'Hans'); }
      catch (_) { /* Keep the normal link working when storage is unavailable. */ }
    });
  });
})();
