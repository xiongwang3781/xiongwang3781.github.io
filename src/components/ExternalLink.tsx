import type {AnchorHTMLAttributes, ReactNode} from 'react';

type ExternalLinkProps = AnchorHTMLAttributes<HTMLAnchorElement> & {
  children: ReactNode;
};

export const ExternalLink = ({children, ...props}: ExternalLinkProps) => (
  <a {...props} target="_blank" rel="noreferrer">
    {children}
    <span className="external-mark" aria-hidden="true">↗</span>
  </a>
);
