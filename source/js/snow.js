/**
 * 下雪特效
 * 浅色模式：淡蓝雪花（白底上可见）；深色模式：白雪花
 * pointer-events: none 不影响页面交互
 */
(function () {
  if (window.matchMedia('(prefers-reduced-motion: reduce)').matches) return;

  var canvas = document.createElement('canvas');
  canvas.id = 'snow-canvas';
  canvas.style.cssText = 'position:fixed;top:0;left:0;width:100%;height:100%;pointer-events:none;z-index:9999;';
  document.body.appendChild(canvas);
  var ctx = canvas.getContext('2d');

  var w, h;
  function resize() {
    w = canvas.width = window.innerWidth;
    h = canvas.height = window.innerHeight;
  }
  resize();
  window.addEventListener('resize', resize);

  var COUNT = window.innerWidth < 768 ? 40 : 90;
  var flakes = [];
  for (var i = 0; i < COUNT; i++) {
    flakes.push({
      x: Math.random() * window.innerWidth,
      y: Math.random() * window.innerHeight,
      r: Math.random() * 2.5 + 1,
      speed: Math.random() * 1 + 0.5,
      swing: Math.random() * 1.5,
      swingStep: Math.random() * Math.PI * 2,
      opacity: Math.random() * 0.5 + 0.4
    });
  }

  function snowColor() {
    return document.documentElement.getAttribute('data-theme') === 'dark'
      ? '255,255,255'
      : '148,184,243';
  }

  function draw() {
    ctx.clearRect(0, 0, w, h);
    var color = snowColor();
    for (var i = 0; i < flakes.length; i++) {
      var f = flakes[i];
      f.swingStep += 0.01;
      f.y += f.speed;
      f.x += Math.sin(f.swingStep) * f.swing * 0.3;
      if (f.y > h + 5) { f.y = -5; f.x = Math.random() * w; }
      if (f.x > w + 5) { f.x = -5; } else if (f.x < -5) { f.x = w + 5; }
      ctx.beginPath();
      ctx.arc(f.x, f.y, f.r, 0, Math.PI * 2);
      ctx.fillStyle = 'rgba(' + color + ',' + f.opacity + ')';
      ctx.fill();
    }
    requestAnimationFrame(draw);
  }
  draw();
})();
