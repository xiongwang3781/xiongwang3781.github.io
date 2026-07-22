# Xiong Wang Academic Site — Brand Specification

## Positioning

An editorial research portfolio for international geoscience peers and potential collaborators. The visual language combines Pentagram-like typographic confidence with Tufte-like information discipline: strong hierarchy, quiet rules, evidence before decoration.

## Design tokens

- Paper: `#FBFAF6`
- Ink: `#1B1B1A`
- Cobalt: `#1E3FFF`
- Data rust: `#A6300E`
- Rule: `#D8D2C2`
- Display: Archivo, self-hosted
- Reading: Source Serif 4, self-hosted
- Grid: 12 columns desktop, 6 tablet, 4 mobile
- Spacing: 8px baseline

## Rules

- Square corners, no drop shadows, no decorative gradients.
- Use cobalt for navigation, actions, and structural emphasis.
- Use rust only for data labels, dates, and research notation.
- Keep scientific prose in a narrow serif measure.
- External links carry a visible north-east arrow.
- Motion is optional, silent, and subordinate to semantic HTML.

## Motion system

- One silent 8-second, 240-frame loop at 30fps.
- Slow photographic scale drift; measured rule reveals; one restrained spring for labels.
- Motion vocabulary: mask, baseline, scale, and opacity. No bounce, glow, or decorative spectacle.
- Static poster is the first paint. The Player is requested during browser idle time.
- Reduced-motion, Save-Data, and low-resource devices retain the static poster and do not load the Player.
