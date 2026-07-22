import type {Profile, Project, Publication, ResearchArea, ResearchOutput} from './types';

export const profile: Profile = {
  name: 'Xiong Wang',
  role: 'Experimental rock physicist',
  institution: 'China University of Geosciences',
  introduction:
    'I study how minerals and rocks deform under high pressure and temperature—and how grain-scale records scale up to the mechanics and seismic character of the lithosphere.',
  email: 'xiong.wang3781@outlook.com',
  cvHref: '/assets/xiong-wang/Xiong-Wang-CV.pdf',
  links: [
    {label: 'Google Scholar', href: 'https://scholar.google.com/citations?hl=en&user=0ugogQ0AAAAJ&view_op=list_works&sortby=pubdate'},
    {label: 'ORCID', href: 'https://orcid.org/0000-0003-3866-880X'},
    {label: 'GitHub', href: 'https://github.com/xiongwang3781'},
  ],
};

export const researchAreas: ResearchArea[] = [
  {
    index: '01',
    title: 'Rock rheology',
    summary: 'Laboratory constraints on the strength and deformation mechanisms of crustal and mantle minerals.',
    keywords: ['HPHT experiments', 'flow laws', 'grain size'],
  },
  {
    index: '02',
    title: 'Mineral microstructures',
    summary: 'Reading fabrics and grain-scale textures as quantitative records of stress, strain, and thermal history.',
    keywords: ['CPO', 'recrystallization', 'texture'],
  },
  {
    index: '03',
    title: 'EBSD & MTEX',
    summary: 'Reproducible crystallographic analysis linking preferred orientation to rheological and seismic properties.',
    keywords: ['EBSD', 'MTEX', 'anisotropy'],
  },
  {
    index: '04',
    title: 'Deep-crust mechanics',
    summary: 'Connecting experimental evidence and natural samples to the mechanical architecture of active orogens.',
    keywords: ['amphibole', 'Tibet', 'lithosphere'],
  },
];

export const publications: Publication[] = [
  {
    year: 2024,
    authors: 'Liu, T., Wang, X., Pu, C. & Jing, Z.C.',
    title: 'Thermoelastic Properties of Seifertite at High Pressures and Temperatures: Implications for Negative Velocity Discontinuities in the D″ Layer',
    journal: 'Geophysical Research Letters 51',
    doi: '10.1029/2024GL112270',
    type: 'Journal article',
  },
  {
    year: 2024,
    authors: 'Liu, Y., Yang, T., Wang, K., Wang, X. & Li, Y.',
    title: 'Influence of Grain Size Evolution on Mantle Plume and LLSVP Dynamics',
    journal: 'Geochemistry, Geophysics, Geosystems 25',
    doi: '10.1029/2024GC011807',
    type: 'Journal article',
  },
  {
    year: 2023,
    authors: 'Wang, X., Zhang, J.F., Tommasi, A., Lopez-Sanchez, M., Jing, Z.C., Shi, F., Liu, W.L. & Barou, F.',
    title: 'Experimental Evidence for a Weak Calcic-Amphibole-Rich Deep Crust in Orogens',
    journal: 'Geophysical Research Letters 50',
    doi: '10.1029/2022GL102320',
    type: 'Journal article',
    featured: true,
  },
  {
    year: 2021,
    authors: 'Wang, X., Zhang, J.F., Tommasi, A., Jing, Z.C. & Yuan, M.S.',
    title: 'Microstructure and Seismic Properties of Amphibole-rich Rocks from the Deep Crust in Southern Tibet',
    journal: 'Tectonophysics 811, 228869',
    doi: '10.1016/j.tecto.2021.228869',
    type: 'Journal article',
  },
  {
    year: 2020,
    authors: 'Li, W.J., Zhang, J.F., Wang, X., Wang, Y., Wu, X. & Hu, Z.C.',
    title: 'Petrofabrics and Seismic Properties of Himalayan Amphibolites: Implications for a Thick Anisotropic Deep Crust Beneath the Tibetan Plateau',
    journal: 'Journal of Geophysical Research: Solid Earth 125',
    doi: '10.1029/2019JB018700',
    type: 'Journal article',
  },
  {
    year: 2019,
    authors: 'Wang, X., Zhang, J.F., Rushmer, T., Adam, J., Turner, S. & Xu, W.',
    title: 'Adakite-Like Potassic Magmatism and Crust–Mantle Interaction in a Postcollisional Setting: An Experimental Study of Melting Beneath the Tibetan Plateau',
    journal: 'Journal of Geophysical Research: Solid Earth 124',
    doi: '10.1029/2019JB018392',
    type: 'Journal article',
  },
];

export const projects: Project[] = [
  {
    funder: 'China Postdoctoral Science Foundation',
    grant: '2024M753023',
    title: 'Rheological properties of plagioclase at high temperature and high pressure',
    role: 'Principal investigator',
    period: 'Jun 2024 — present',
  },
  {
    funder: 'National Natural Science Foundation of China',
    grant: '42302242',
    title: 'Rheological properties of calcic amphibole at high temperature and high pressure',
    role: 'Principal investigator',
    period: 'Jan 2024 — Dec 2026',
  },
  {
    funder: 'China Postdoctoral Science Foundation',
    grant: '2021M701563',
    title: 'Calcic-amphibole-rich deep crust beneath southern Tibet',
    role: 'Principal investigator',
    period: 'Oct 2021 — Aug 2022',
    status: 'Completed',
  },
];

export const researchOutputs: ResearchOutput[] = [
  {year: 2026, type: 'Supplementary material', title: 'Supplementary material for a Geology research contribution', doi: '10.1130/GEOL.S.33006923.v1'},
  {year: 2025, type: 'Dataset', title: 'Open research dataset', doi: '10.6084/m9.figshare.29299376'},
  {year: 2025, type: 'Dataset', title: 'Open research dataset', doi: '10.6084/m9.figshare.28120049'},
  {year: 2022, type: 'Preprint', title: 'Earth and Space Science Open Archive preprint', doi: '10.1002/essoar.10510645.1'},
  {year: 2022, type: 'Dataset', title: 'Open research dataset', doi: '10.6084/m9.figshare.20227329'},
];
