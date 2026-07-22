import {lazy, Suspense, useCallback, useEffect, useMemo, useRef, useState} from 'react';
import type {MouseEvent} from 'react';
import type {PlayerRef} from '@remotion/player';
import type {ResearchScaleHeroProps} from '../remotion/ResearchScaleHero';

const LazyHeroPlayer = lazy(() =>
  import('./HeroPlayer').then((module) => ({default: module.HeroPlayer})),
);

type MotionMode = 'loading' | 'player' | 'reduced' | 'save-data' | 'low-performance';

type NavigatorWithPerformanceHints = Navigator & {
  connection?: {saveData?: boolean};
  deviceMemory?: number;
};

const getFallbackMode = (): Exclude<MotionMode, 'loading' | 'player'> | null => {
  const browserNavigator = navigator as NavigatorWithPerformanceHints;
  if (window.matchMedia('(prefers-reduced-motion: reduce)').matches) return 'reduced';
  if (browserNavigator.connection?.saveData) return 'save-data';
  if (
    (browserNavigator.deviceMemory !== undefined && browserNavigator.deviceMemory <= 2) ||
    navigator.hardwareConcurrency <= 2
  ) return 'low-performance';
  return null;
};

const motionStatusCopy: Record<MotionMode, string> = {
  loading: 'STATIC POSTER / PLAYER LOADING',
  player: 'REMOTION / 8 SECOND LOOP',
  reduced: 'STATIC / REDUCED MOTION',
  'save-data': 'STATIC / SAVE-DATA',
  'low-performance': 'STATIC / LIGHTWEIGHT MODE',
};

const HeroPoster = () => (
  <div className="hero-poster" aria-hidden="true">
    <picture>
      <source media="(max-width: 760px)" srcSet="/assets/xiong-wang/field-mobile-1200.avif" type="image/avif" />
      <source media="(max-width: 760px)" srcSet="/assets/xiong-wang/field-mobile-1200.webp" type="image/webp" />
      <source srcSet="/assets/xiong-wang/field-2400.avif" type="image/avif" />
      <img src="/assets/xiong-wang/field-2400.webp" alt="" />
    </picture>
    <div className="hero-poster-veil" />
    <div className="hero-poster-rules">{Array.from({length: 9}, (_, index) => <i key={index} />)}</div>
    <div className="hero-poster-scale">
      <span>mineral</span><span>grain</span><span>rock</span><span>lithosphere</span>
    </div>
    <div className="hero-poster-caption">FIELD / LAB / CODE</div>
  </div>
);

export const HeroVisual = () => {
  const containerRef = useRef<HTMLDivElement>(null);
  const playerRef = useRef<PlayerRef>(null);
  const userPausedRef = useRef(false);
  const resumeOnReturnRef = useRef(false);
  const [motionMode, setMotionMode] = useState<MotionMode>('loading');
  const [shouldLoadPlayer, setShouldLoadPlayer] = useState(false);
  const [playerReady, setPlayerReady] = useState(false);
  const [isPlaying, setIsPlaying] = useState(false);
  const inputProps = useMemo<ResearchScaleHeroProps>(() => ({
    imageSrc: '/assets/xiong-wang/field-2400.avif',
  }), []);

  const handleReady = useCallback(() => {
    setPlayerReady(true);
    setMotionMode('player');
    setIsPlaying(playerRef.current?.isPlaying() ?? true);
  }, []);

  useEffect(() => {
    const fallbackMode = getFallbackMode();
    if (fallbackMode) {
      setMotionMode(fallbackMode);
      return;
    }

    const loadPlayer = () => setShouldLoadPlayer(true);
    if ('requestIdleCallback' in window) {
      const idleId = window.requestIdleCallback(loadPlayer, {timeout: 1200});
      return () => window.cancelIdleCallback(idleId);
    }
    const timeoutId = globalThis.setTimeout(loadPlayer, 450);
    return () => globalThis.clearTimeout(timeoutId);
  }, []);

  useEffect(() => {
    const mediaQuery = window.matchMedia('(prefers-reduced-motion: reduce)');
    const handlePreferenceChange = (event: MediaQueryListEvent) => {
      if (!event.matches) return;
      playerRef.current?.pause();
      setIsPlaying(false);
      setShouldLoadPlayer(false);
      setPlayerReady(false);
      setMotionMode('reduced');
    };
    mediaQuery.addEventListener('change', handlePreferenceChange);
    return () => mediaQuery.removeEventListener('change', handlePreferenceChange);
  }, []);

  useEffect(() => {
    if (!playerReady || !playerRef.current) return;
    const player = playerRef.current;
    const handlePlay = () => setIsPlaying(true);
    const handlePause = () => setIsPlaying(false);
    player.addEventListener('play', handlePlay);
    player.addEventListener('pause', handlePause);
    return () => {
      player.removeEventListener('play', handlePlay);
      player.removeEventListener('pause', handlePause);
    };
  }, [playerReady]);

  useEffect(() => {
    if (!playerReady || !containerRef.current) return;
    const observer = new IntersectionObserver(([entry]) => {
      const player = playerRef.current;
      if (!player) return;
      if (entry.isIntersecting && entry.intersectionRatio >= 0.18) {
        if (resumeOnReturnRef.current && !userPausedRef.current) player.play();
        resumeOnReturnRef.current = false;
        return;
      }
      if (player.isPlaying() && !userPausedRef.current) resumeOnReturnRef.current = true;
      player.pause();
    }, {threshold: [0, 0.18, 0.5]});
    observer.observe(containerRef.current);
    return () => observer.disconnect();
  }, [playerReady]);

  const toggleMotion = (event: MouseEvent<HTMLButtonElement>) => {
    const player = playerRef.current;
    if (!player) return;
    if (player.isPlaying()) {
      userPausedRef.current = true;
      player.pause();
    } else {
      userPausedRef.current = false;
      player.play(event);
    }
  };

  return (
    <div
      ref={containerRef}
      className="hero-visual"
      data-motion-mode={motionMode}
      data-motion-playing={isPlaying}
    >
      <HeroPoster />
      {shouldLoadPlayer && (
        <Suspense fallback={null}>
          <div className="hero-player-layer" data-ready={playerReady} aria-hidden="true" inert>
            <div>
              <LazyHeroPlayer playerRef={playerRef} inputProps={inputProps} onReady={handleReady} />
            </div>
          </div>
        </Suspense>
      )}
      {motionMode === 'player' && (
        <button className="motion-control" type="button" onClick={toggleMotion}>
          <span className="motion-control-mark" aria-hidden="true">{isPlaying ? 'Ⅱ' : '▶'}</span>
          {isPlaying ? 'Pause motion' : 'Play motion'}
        </button>
      )}
      <span className="motion-status">{motionStatusCopy[motionMode]}</span>
    </div>
  );
};
