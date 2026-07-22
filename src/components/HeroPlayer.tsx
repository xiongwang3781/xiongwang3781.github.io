import {useEffect} from 'react';
import type {RefObject} from 'react';
import {Player} from '@remotion/player';
import type {PlayerRef} from '@remotion/player';
import {ResearchScaleHero} from '../remotion/ResearchScaleHero';
import type {ResearchScaleHeroProps} from '../remotion/ResearchScaleHero';

type HeroPlayerProps = {
  playerRef: RefObject<PlayerRef | null>;
  inputProps: ResearchScaleHeroProps;
  onReady: () => void;
};

export const HeroPlayer = ({playerRef, inputProps, onReady}: HeroPlayerProps) => {
  useEffect(() => {
    const animationFrame = requestAnimationFrame(onReady);
    return () => cancelAnimationFrame(animationFrame);
  }, [onReady]);

  return (
    <Player
      ref={playerRef}
      component={ResearchScaleHero}
      inputProps={inputProps}
      durationInFrames={240}
      compositionWidth={1920}
      compositionHeight={1080}
      fps={30}
      autoPlay
      loop
      initiallyMuted
      controls={false}
      showVolumeControls={false}
      clickToPlay={false}
      spaceKeyToPlayOrPause={false}
      acknowledgeRemotionLicense
      style={{width: '100%'}}
    />
  );
};
