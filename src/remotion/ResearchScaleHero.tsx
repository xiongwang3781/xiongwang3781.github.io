import {
  AbsoluteFill,
  Easing,
  Img,
  Sequence,
  interpolate,
  spring,
  staticFile,
  useCurrentFrame,
  useVideoConfig,
} from 'remotion';

export type ResearchScaleHeroProps = {
  imageSrc?: string;
};

const scaleLabels = ['mineral', 'grain', 'rock', 'lithosphere'];

const MeasurementGrid = () => {
  const frame = useCurrentFrame();
  const lineProgress = interpolate(frame, [0, 42], [0, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
    easing: Easing.out(Easing.cubic),
  });

  return (
    <div className="hero-composition-rules">
      {Array.from({length: 9}, (_, index) => {
        const stagger = Math.max(0, Math.min(1, lineProgress * 1.45 - index * 0.055));
        return (
          <i
            key={index}
            style={{
              opacity: stagger * 0.68,
              transform: `scaleY(${0.12 + stagger * 0.88})`,
              transformOrigin: index % 2 === 0 ? 'bottom' : 'top',
            }}
          />
        );
      })}
    </div>
  );
};

const ResearchScale = () => {
  const frame = useCurrentFrame();
  const {fps} = useVideoConfig();
  const settle = spring({
    fps,
    frame,
    config: {damping: 22, stiffness: 90, mass: 0.8},
    durationInFrames: 42,
  });

  return (
    <div
      className="hero-composition-scale"
      style={{
        opacity: settle,
        transform: `translateY(${interpolate(settle, [0, 1], [34, 0])}px)`,
      }}
    >
      {scaleLabels.map((label, index) => {
        const labelOpacity = interpolate(frame, [index * 8, index * 8 + 16], [0, 1], {
          extrapolateLeft: 'clamp',
          extrapolateRight: 'clamp',
        });
        return <span key={label} style={{opacity: labelOpacity}}>{label}</span>;
      })}
    </div>
  );
};

const ScaleStatement = () => {
  const frame = useCurrentFrame();
  const {fps} = useVideoConfig();
  const entry = spring({
    fps,
    frame,
    config: {damping: 20, stiffness: 82, mass: 0.9},
    durationInFrames: 48,
  });
  const sweep = interpolate(frame, [20, 92], [0, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
    easing: Easing.inOut(Easing.cubic),
  });

  return (
    <div className="hero-composition-statement">
      <div
        className="hero-composition-statement-rule"
        style={{transform: `scaleX(${sweep})`}}
      />
      <div
        className="hero-composition-statement-copy"
        style={{
          opacity: entry,
          transform: `translateY(${interpolate(entry, [0, 1], [52, 0])}px)`,
        }}
      >
        <span>FORCE</span>
        <span className="statement-arrow">→</span>
        <span>FABRIC</span>
      </div>
    </div>
  );
};

const FieldMarker = () => {
  const frame = useCurrentFrame();
  const reveal = interpolate(frame, [0, 24, 84, 110], [0, 1, 1, 0], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });
  const x = interpolate(frame, [0, 110], [14, 72], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
    easing: Easing.inOut(Easing.cubic),
  });

  return (
    <div
      className="hero-composition-marker"
      style={{left: `${x}%`, opacity: reveal}}
    >
      <i />
      <span>sample / 04</span>
    </div>
  );
};

export const ResearchScaleHero = ({
  imageSrc = staticFile('assets/xiong-wang/field-2400.avif'),
}: ResearchScaleHeroProps) => {
  const frame = useCurrentFrame();
  const cycleFade = interpolate(frame, [0, 18, 212, 239], [0, 1, 1, 0], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
    easing: Easing.inOut(Easing.cubic),
  });
  const imageScale = interpolate(frame, [0, 120, 239], [1.018, 1.064, 1.018], {
    easing: Easing.inOut(Easing.sin),
  });
  const imageX = interpolate(frame, [0, 120, 239], [-0.7, 0.7, -0.7], {
    easing: Easing.inOut(Easing.sin),
  });
  const maskReveal = interpolate(frame, [6, 58, 184, 230], [18, 96, 96, 18], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
    easing: Easing.inOut(Easing.cubic),
  });

  return (
    <AbsoluteFill className="hero-composition">
      <Img
        className="hero-composition-image hero-composition-image-base"
        src={imageSrc}
        style={{transform: `translateX(${imageX}%) scale(${imageScale})`}}
      />
      <Img
        className="hero-composition-image hero-composition-image-focus"
        src={imageSrc}
        style={{
          clipPath: `inset(0 ${100 - maskReveal}% 0 0)`,
          transform: `translateX(${imageX}%) scale(${imageScale})`,
        }}
      />
      <AbsoluteFill className="hero-composition-veil" style={{opacity: cycleFade}} />
      <Sequence from={12} durationInFrames={198} premountFor={12}>
        <MeasurementGrid />
      </Sequence>
      <Sequence from={28} durationInFrames={182} premountFor={12}>
        <ResearchScale />
      </Sequence>
      <Sequence from={64} durationInFrames={132} premountFor={10}>
        <ScaleStatement />
      </Sequence>
      <Sequence from={36} durationInFrames={126} premountFor={8}>
        <FieldMarker />
      </Sequence>
      <div className="hero-composition-caption" style={{opacity: cycleFade}}>
        FIELD / LAB / CODE
      </div>
      <div className="hero-composition-frame-count" style={{opacity: cycleFade}}>
        {String(frame + 1).padStart(3, '0')} / 240
      </div>
    </AbsoluteFill>
  );
};
