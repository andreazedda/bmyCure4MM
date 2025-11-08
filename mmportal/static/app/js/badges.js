(function (root) {
  function toNumber(value) {
    var num = parseFloat(value);
    return Number.isFinite(num) ? num : NaN;
  }

  function computeDoseStatus(value, min, max) {
    var v = toNumber(value);
    var lo = toNumber(min);
    var hi = toNumber(max);
    if (!Number.isFinite(v) || !Number.isFinite(lo) || !Number.isFinite(hi) || lo >= hi) {
      return 'unknown';
    }
    if (v < lo || v > hi) {
      return 'danger';
    }
    var warnBand = Math.max((hi - lo) * 0.05, 0.5);
    if (v - lo < warnBand || hi - v < warnBand) {
      return 'warn';
    }
    return 'safe';
  }

  function computeInteractionStatus(value) {
    var v = toNumber(value);
    if (!Number.isFinite(v)) {
      return 'unknown';
    }
    if (v > 0.2) {
      return 'danger';
    }
    if (v >= 0.1) {
      return 'warn';
    }
    return 'safe';
  }

  function debounce(fn, ms) {
    ms = ms || 150;
    var timeout;
    return function () {
      var args = arguments;
      var context = this;
      clearTimeout(timeout);
      timeout = setTimeout(function () {
        fn.apply(context, args);
      }, ms);
    };
  }

  var exportsObj = {
    computeDoseStatus: computeDoseStatus,
    computeInteractionStatus: computeInteractionStatus,
    debounce: debounce,
  };
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = exportsObj;
  } else {
    root.BadgeUtils = exportsObj;
  }
})(typeof self !== 'undefined' ? self : globalThis);
