(() => {
  const pokeball = document.querySelector('.floating-pokeball');
  if (!pokeball) return;

  const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)');
  let animationFrame = 0;

  const clearRotation = () => {
    pokeball.style.removeProperty('--ball-turn');
  };

  const rotatePokeball = () => {
    animationFrame = 0;

    if (reducedMotion.matches) {
      clearRotation();
      return;
    }

    const scrollY = Math.max(window.scrollY, 0);
    pokeball.style.setProperty('--ball-turn', `${scrollY * 0.35}deg`);
  };

  const scheduleRotation = () => {
    if (reducedMotion.matches || animationFrame) return;
    animationFrame = window.requestAnimationFrame(rotatePokeball);
  };

  const handleMotionPreference = () => {
    if (animationFrame) {
      window.cancelAnimationFrame(animationFrame);
      animationFrame = 0;
    }

    clearRotation();
    scheduleRotation();
  };

  window.addEventListener('scroll', scheduleRotation, { passive: true });
  window.addEventListener('pageshow', scheduleRotation);
  reducedMotion.addEventListener('change', handleMotionPreference);

  scheduleRotation();
})();
