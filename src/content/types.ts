export type ExternalLink = {
  label: string;
  href: string;
};

export type Profile = {
  name: string;
  role: string;
  institution: string;
  introduction: string;
  email: string;
  cvHref: string;
  links: ExternalLink[];
};

export type ResearchArea = {
  index: string;
  title: string;
  summary: string;
  keywords: string[];
};

export type Project = {
  funder: string;
  grant: string;
  title: string;
  role: string;
  period: string;
  status?: string;
};

export type Publication = {
  year: number;
  authors: string;
  title: string;
  journal: string;
  doi: string;
  type: 'Journal article';
  featured?: boolean;
};

export type ResearchOutput = {
  year: number;
  type: 'Dataset' | 'Preprint' | 'Supplementary material';
  title: string;
  doi: string;
};
