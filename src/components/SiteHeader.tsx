import {useState} from 'react';

const navigation = [
  ['Research', '#research'],
  ['Publications', '#publications'],
  ['Projects', '#projects'],
  ['Code', '#code'],
  ['Profile', '#profile'],
] as const;

export const SiteHeader = () => {
  const [open, setOpen] = useState(false);

  return (
    <header className="site-header">
      <a className="wordmark" href="#top" aria-label="Xiong Wang, home">
        <span className="wordmark-monogram" aria-hidden="true">XW</span>
        <span className="wordmark-name">Xiong Wang</span>
      </a>
      <button
        className="menu-toggle"
        type="button"
        aria-expanded={open}
        aria-controls="primary-navigation"
        onClick={() => setOpen((current) => !current)}
      >
        <span>{open ? 'Close' : 'Menu'}</span>
        <span className="menu-glyph" aria-hidden="true">{open ? '×' : '＋'}</span>
      </button>
      <nav id="primary-navigation" className="site-nav" data-open={open} aria-label="Primary navigation">
        {navigation.map(([label, href]) => (
          <a key={href} href={href} onClick={() => setOpen(false)}>{label}</a>
        ))}
      </nav>
    </header>
  );
};
