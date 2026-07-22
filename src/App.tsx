import {ExternalLink} from './components/ExternalLink';
import {HeroVisual} from './components/HeroVisual';
import {SiteHeader} from './components/SiteHeader';
import {profile, projects, publications, researchAreas, researchOutputs} from './content/siteContent';

const doiHref = (doi: string) => `https://doi.org/${doi}`;

export const App = () => {
  const featured = publications.find((publication) => publication.featured)!;
  const remainingPublications = publications.filter((publication) => !publication.featured);

  return (
    <>
      <a className="skip-link" href="#main-content">Skip to content</a>
      <SiteHeader />
      <main id="main-content">
        <section className="hero section-shell" id="top" aria-labelledby="hero-title">
          <div className="hero-copy">
            <p className="kicker">Experimental rock physics / Wuhan, China</p>
            <h1 id="hero-title">How rocks carry the memory of force.</h1>
            <p className="hero-deck">{profile.introduction}</p>
            <div className="hero-actions">
              <a className="action action-primary" href={`mailto:${profile.email}`}>Start a conversation</a>
              <a className="action action-secondary" href="#research">Explore research</a>
            </div>
          </div>
          <HeroVisual />
          <aside className="hero-note" aria-label="Research approach">
            <span className="note-index">01—04</span>
            <p>Experiments, microstructures, crystallography, and computation across scales.</p>
          </aside>
        </section>

        <section className="section section-shell" id="research" aria-labelledby="research-title">
          <header className="section-heading">
            <p className="kicker">Research</p>
            <h2 id="research-title">From lattice orientation to lithosphere mechanics.</h2>
            <p>Four connected lines of inquiry organize the work—not as isolated topics, but as measurements at different scales.</p>
          </header>
          <div className="research-grid">
            {researchAreas.map((area) => (
              <article className="research-area" key={area.index}>
                <span className="area-index">{area.index}</span>
                <h3>{area.title}</h3>
                <p>{area.summary}</p>
                <ul className="keyword-list" aria-label={`${area.title} keywords`}>
                  {area.keywords.map((keyword) => <li key={keyword}>{keyword}</li>)}
                </ul>
              </article>
            ))}
          </div>
        </section>

        <section className="section publication-section" id="publications" aria-labelledby="publications-title">
          <div className="section-shell">
            <header className="section-heading split-heading">
              <div>
                <p className="kicker">Selected publications</p>
                <h2 id="publications-title">Evidence, tested under pressure.</h2>
              </div>
              <ExternalLink href={profile.links[0].href}>Complete record on Google Scholar</ExternalLink>
            </header>
            <article className="featured-publication">
              <div className="featured-label">Featured / {featured.year}</div>
              <div>
                <h3>{featured.title}</h3>
                <p className="publication-authors">{featured.authors}</p>
                <p className="publication-journal">{featured.journal}</p>
                <ExternalLink href={doiHref(featured.doi)}>DOI {featured.doi}</ExternalLink>
              </div>
            </article>
            <ol className="publication-list">
              {remainingPublications.map((publication) => (
                <li key={publication.doi}>
                  <time>{publication.year}</time>
                  <div>
                    <h3><ExternalLink href={doiHref(publication.doi)}>{publication.title}</ExternalLink></h3>
                    <p className="publication-authors">{publication.authors}</p>
                    <p className="publication-journal">{publication.journal}</p>
                  </div>
                  <span className="publication-type">Article</span>
                </li>
              ))}
            </ol>
          </div>
        </section>

        <section className="section section-shell" id="projects" aria-labelledby="projects-title">
          <header className="section-heading">
            <p className="kicker">Funded projects</p>
            <h2 id="projects-title">Questions with institutional support.</h2>
          </header>
          <div className="project-list">
            {projects.map((project, index) => (
              <article className="project" key={project.grant}>
                <span className="project-index">P{String(index + 1).padStart(2, '0')}</span>
                <div className="project-main">
                  <p className="project-funder">{project.funder}</p>
                  <h3>{project.title}</h3>
                </div>
                <dl className="project-meta">
                  <div><dt>Grant</dt><dd>{project.grant}</dd></div>
                  <div><dt>Role</dt><dd>{project.role}</dd></div>
                  <div><dt>Period</dt><dd>{project.period}</dd></div>
                  {project.status && <div><dt>Status</dt><dd>{project.status}</dd></div>}
                </dl>
              </article>
            ))}
          </div>
        </section>

        <section className="section code-section" id="code" aria-labelledby="code-title">
          <div className="section-shell code-grid">
            <header>
              <p className="kicker kicker-light">AIMER / Research code</p>
              <h2 id="code-title">Methods should travel with the result.</h2>
            </header>
            <div className="code-copy">
              <p>AIMER is a developing research-code direction connecting domain knowledge, microstructure analysis, and AI-assisted programming. The current repository includes working MATLAB scripts rather than a packaged software product.</p>
              <ul className="file-list">
                <li><code>Deformation_mechanism_map_ol.m</code><span>Dry-olivine mechanism maps</span></li>
                <li><code>Griggs_data_reduction_cug_v2.m</code><span>Experimental data reduction</span></li>
                <li><code>Diffusion_time_fitting.m</code><span>Diffusion-time fitting</span></li>
                <li><code>elastic_vein_model_rubin1995.m</code><span>Elastic vein model</span></li>
              </ul>
              <ExternalLink className="action action-on-dark" href="https://github.com/xiongwang3781/xiongwang3781.github.io/tree/main/AIMER">Inspect the MATLAB source</ExternalLink>
            </div>
          </div>
        </section>

        <section className="section section-shell outputs-section" aria-labelledby="outputs-title">
          <header className="section-heading split-heading">
            <div>
              <p className="kicker">Open research outputs</p>
              <h2 id="outputs-title">Datasets and materials, counted separately.</h2>
            </div>
            <p>These records are intentionally not mixed with the peer-reviewed publication list.</p>
          </header>
          <ol className="output-list">
            {researchOutputs.map((output) => (
              <li key={output.doi}>
                <time>{output.year}</time>
                <span>{output.type}</span>
                <ExternalLink href={doiHref(output.doi)}>{output.title}</ExternalLink>
              </li>
            ))}
          </ol>
        </section>

        <section className="section profile-section" id="profile" aria-labelledby="profile-title">
          <div className="section-shell profile-grid">
            <figure className="portrait-frame">
              <picture>
                <source srcSet="/assets/xiong-wang/profile-640.avif" type="image/avif" />
                <img src="/assets/xiong-wang/profile-640.webp" alt="Portrait of Xiong Wang" width="358" height="441" />
              </picture>
              <figcaption>Experimental rock physicist<br />Wuhan, China</figcaption>
            </figure>
            <div className="profile-copy">
              <p className="kicker">Profile & contact</p>
              <h2 id="profile-title">Xiong Wang</h2>
              <p className="profile-role">{profile.role}<br />{profile.institution}</p>
              <p>I welcome conversations about experimental rock mechanics, crystallographic analysis, reproducible research workflows, and collaborations that connect laboratory evidence with geodynamics.</p>
              <div className="profile-actions">
                <a className="action action-primary" href={`mailto:${profile.email}`}>{profile.email}</a>
                <a className="action action-secondary" href={profile.cvHref} download>Download CV</a>
              </div>
              <div className="profile-links" aria-label="Academic profiles">
                {profile.links.map((link) => <ExternalLink href={link.href} key={link.label}>{link.label}</ExternalLink>)}
              </div>
            </div>
          </div>
        </section>
      </main>
      <footer className="site-footer section-shell">
        <p>© 2026 Xiong Wang</p>
        <p>Built as an open, static research record.</p>
        <a href="#top">Back to top ↑</a>
      </footer>
    </>
  );
};
