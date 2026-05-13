'use strict';

// ─── Constants ────────────────────────────────────────────────────────────────

const FREQ_MIN   = 20;
const FREQ_MAX   = 20000;
const N_POINTS   = 400;

// Foobar2000 標準 EQ (.feq) の18バンド周波数
const FEQ_FREQS = [
  55, 77, 110, 156, 220, 311, 440, 622,
  880, 1200, 1800, 2500, 3500, 5000,
  7000, 10000, 14000, 20000,
];

// XGEQ エクスポート用 ISO 1/3オクターブ バンド（31バンド）
const XGEQ_EXPORT_FREQS = [
  20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160,
  200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600,
  2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000,
];

// Pre-computed log-spaced frequency array for chart rendering
const freqPoints = logspace(FREQ_MIN, FREQ_MAX, N_POINTS);

// ─── Math Utilities ───────────────────────────────────────────────────────────

function logspace(start, end, n) {
  const a = Math.log10(start);
  const b = Math.log10(end);
  return Array.from({ length: n }, (_, i) => Math.pow(10, a + (b - a) * i / (n - 1)));
}

// ─── Natural Cubic Spline (in log-frequency space) ───────────────────────────

class CubicSpline {
  constructor(xs, ys) {
    const n = xs.length;
    this.xs = xs;
    this.ys = ys;
    this.n  = n;

    if (n < 3) { this.b = this.c = this.d = null; return; }

    const h = Array.from({ length: n - 1 }, (_, i) => xs[i + 1] - xs[i]);

    // Tridiagonal system for second derivatives (natural BCs: c[0] = c[n-1] = 0)
    const alpha = new Array(n).fill(0);
    for (let i = 1; i < n - 1; i++) {
      alpha[i] = (3 / h[i]) * (ys[i + 1] - ys[i]) - (3 / h[i - 1]) * (ys[i] - ys[i - 1]);
    }

    const l  = new Array(n).fill(1);
    const mu = new Array(n).fill(0);
    const z  = new Array(n).fill(0);

    for (let i = 1; i < n - 1; i++) {
      l[i]  = 2 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1];
      mu[i] = h[i] / l[i];
      z[i]  = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    const c = new Array(n).fill(0);
    const b = new Array(n - 1);
    const d = new Array(n - 1);

    for (let j = n - 2; j >= 0; j--) {
      c[j] = z[j] - mu[j] * c[j + 1];
      b[j] = (ys[j + 1] - ys[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
      d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    this.b = b; this.c = c; this.d = d;
  }

  evaluate(x) {
    const { xs, ys, n, b, c, d } = this;
    if (n === 0) return 0;
    if (n === 1) return ys[0];
    if (x <= xs[0])     return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];

    // Binary search for segment
    let lo = 0, hi = n - 2;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (xs[mid] <= x) lo = mid; else hi = mid - 1;
    }
    const i  = lo;
    const dx = x - xs[i];

    if (!b) {
      // Fallback: linear interpolation for n === 2
      return ys[i] + (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]) * dx;
    }

    return ys[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
  }
}

// ─── EQ Interpolation ─────────────────────────────────────────────────────────

function buildSpline(bands) {
  const sorted = bands
    .filter(b => b.freq > 0 && isFinite(b.freq) && isFinite(b.gain))
    .sort((a, b) => a.freq - b.freq);

  // Deduplicate by frequency (keep last)
  const deduped = [];
  for (const band of sorted) {
    if (deduped.length && deduped[deduped.length - 1].freq === band.freq) {
      deduped[deduped.length - 1] = band;
    } else {
      deduped.push(band);
    }
  }

  return deduped;
}

function interpolateEQ(bands) {
  const pts = buildSpline(bands);

  if (pts.length === 0) return freqPoints.map(() => 0);
  if (pts.length === 1) return freqPoints.map(() => pts[0].gain);

  const lxs = pts.map(p => Math.log10(p.freq));
  const ys  = pts.map(p => p.gain);

  if (pts.length === 2) {
    return freqPoints.map(f => {
      const lx = Math.log10(f);
      if (lx <= lxs[0]) return ys[0];
      if (lx >= lxs[1]) return ys[1];
      return ys[0] + (ys[1] - ys[0]) * (lx - lxs[0]) / (lxs[1] - lxs[0]);
    });
  }

  const spline = new CubicSpline(lxs, ys);
  return freqPoints.map(f => spline.evaluate(Math.log10(f)));
}

function gainAtFreq(bands, freq) {
  const pts = buildSpline(bands);
  if (pts.length === 0) return 0;
  if (pts.length === 1) return pts[0].gain;

  const lx  = Math.log10(Math.max(FREQ_MIN, Math.min(FREQ_MAX, freq)));
  const lxs = pts.map(p => Math.log10(p.freq));
  const ys  = pts.map(p => p.gain);

  if (pts.length === 2) {
    if (lx <= lxs[0]) return ys[0];
    if (lx >= lxs[1]) return ys[1];
    return ys[0] + (ys[1] - ys[0]) * (lx - lxs[0]) / (lxs[1] - lxs[0]);
  }

  return new CubicSpline(lxs, ys).evaluate(lx);
}

// ─── State ────────────────────────────────────────────────────────────────────

const state = {
  eq1: [
    { freq: 60,    gain:  3.0 },
    { freq: 250,   gain: -2.0 },
    { freq: 1000,  gain:  0.0 },
    { freq: 4000,  gain: -1.0 },
    { freq: 12000, gain:  3.0 },
  ],
  eq2: [
    { freq: 80,   gain:  5.0 },
    { freq: 500,  gain:  1.0 },
    { freq: 2000, gain: -2.0 },
    { freq: 8000, gain:  2.0 },
  ],
};

let decimalPlaces = 2;

// ─── Chart ────────────────────────────────────────────────────────────────────

let chart;

function initChart() {
  const ctx = document.getElementById('eq-chart').getContext('2d');

  chart = new Chart(ctx, {
    type: 'line',
    data: {
      datasets: [
        {
          label: 'EQ 1',
          data: [],
          borderColor: '#818cf8',
          borderWidth: 2,
          borderDash: [6, 4],
          pointRadius: 0,
          tension: 0,
        },
        {
          label: 'EQ 2',
          data: [],
          borderColor: '#f472b6',
          borderWidth: 2,
          borderDash: [6, 4],
          pointRadius: 0,
          tension: 0,
        },
        {
          label: 'Mixed',
          data: [],
          borderColor: '#34d399',
          borderWidth: 3,
          borderDash: [],
          pointRadius: 0,
          tension: 0,
          fill: {
            target: 'origin',
            above: 'rgba(52, 211, 153, 0.05)',
            below: 'rgba(52, 211, 153, 0.05)',
          },
        },
      ],
    },
    options: {
      animation: { duration: 120 },
      responsive: true,
      maintainAspectRatio: false,
      interaction: { mode: 'index', intersect: false },
      onClick(event, _elements, ch) {
        // Chart.js v4: event is ChartEvent, .x/.y are already canvas-relative
        if (event.x == null) return;
        const freq = ch.scales.x.getValueForPixel(event.x);
        if (freq && isFinite(freq) && freq > 0) addQueryFromChart(freq);
      },
      plugins: {
        legend: { display: false },
        tooltip: {
          backgroundColor: 'rgba(8, 12, 28, 0.92)',
          borderColor: 'rgba(255,255,255,0.1)',
          borderWidth: 1,
          titleColor: '#94a3b8',
          bodyColor: '#f1f5f9',
          padding: 10,
          callbacks: {
            title([item]) {
              const f = item.parsed.x;
              return f >= 1000
                ? `${(f / 1000).toPrecision(3).replace(/\.?0+$/, '')} kHz`
                : `${Math.round(f)} Hz`;
            },
            label(item) {
              const sign = item.parsed.y >= 0 ? '+' : '';
              return `  ${item.dataset.label}: ${sign}${item.parsed.y.toFixed(2)} dB`;
            },
          },
        },
      },
      scales: {
        x: {
          type: 'logarithmic',
          min: FREQ_MIN,
          max: FREQ_MAX,
          grid: {
            color: 'rgba(255,255,255,0.04)',
          },
          border: {
            color: 'rgba(255,255,255,0.08)',
          },
          ticks: {
            color: '#475569',
            maxRotation: 0,
            autoSkip: false,
            callback(val) {
              const marks = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000];
              const hit = marks.find(m => Math.abs(m - val) / m < 0.015);
              if (!hit) return null;
              return hit >= 1000 ? hit / 1000 + 'k' : String(hit);
            },
          },
        },
        y: {
          min: -24,
          max: 24,
          grid: {
            color: (ctx) => ctx.tick?.value === 0
              ? 'rgba(255,255,255,0.14)'
              : 'rgba(255,255,255,0.04)',
          },
          border: {
            color: 'rgba(255,255,255,0.08)',
            dash: [4, 4],
          },
          ticks: {
            color: '#475569',
            stepSize: 6,
            callback: v => (v > 0 ? '+' : '') + v + ' dB',
          },
        },
      },
    },
  });
}

function updateChart() {
  if (!chart) return;
  const y1  = interpolateEQ(state.eq1);
  const y2  = interpolateEQ(state.eq2);
  const ym  = y1.map((v, i) => (v + y2[i]) / 2);
  const toXY = vals => freqPoints.map((f, i) => ({ x: f, y: vals[i] }));

  chart.data.datasets[0].data = toXY(y1);
  chart.data.datasets[1].data = toXY(y2);
  chart.data.datasets[2].data = toXY(ym);
  chart.update('none');
}

// ─── UI: Band Rows ────────────────────────────────────────────────────────────

function renderBands(eqKey) {
  const bands     = state[eqKey];
  const container = document.getElementById(`bands-${eqKey}`);
  container.innerHTML = '';

  bands.forEach((band, idx) => {
    const row = document.createElement('div');
    row.className = 'band-row';

    row.innerHTML = `
      <div class="input-group">
        <input
          class="band-input"
          type="number"
          value="${band.freq}"
          min="1" max="96000" step="1"
          data-eq="${eqKey}" data-idx="${idx}" data-field="freq"
          aria-label="周波数 ${idx + 1}"
          placeholder="1000"
        >
        <span class="input-unit">Hz</span>
      </div>
      <div class="input-group">
        <input
          class="band-input"
          type="number"
          value="${band.gain}"
          min="-40" max="40" step="0.1"
          data-eq="${eqKey}" data-idx="${idx}" data-field="gain"
          aria-label="ゲイン ${idx + 1}"
          placeholder="0"
        >
        <span class="input-unit">dB</span>
      </div>
      <button
        class="btn-remove"
        data-eq="${eqKey}" data-idx="${idx}"
        aria-label="バンド ${idx + 1} を削除"
        title="削除"
      >
        <svg width="12" height="12" viewBox="0 0 12 12" fill="none">
          <path d="M1 1l10 10M11 1L1 11" stroke="currentColor" stroke-width="2" stroke-linecap="round"/>
        </svg>
      </button>
    `;

    container.appendChild(row);
  });
}

function addBand(eqKey) {
  state[eqKey].push({ freq: 1000, gain: 0 });
  renderBands(eqKey);
  updateChart();
}

function removeBand(eqKey, idx) {
  state[eqKey].splice(idx, 1);
  renderBands(eqKey);
  updateChart();
}

function updateBand(eqKey, idx, field, raw) {
  const val = parseFloat(raw);
  if (!isFinite(val)) return;
  state[eqKey][idx][field] = val;
  updateChart();
}

function clearEQ(eqKey) {
  state[eqKey] = [];
  renderBands(eqKey);
  updateChart();
}

// ─── Query Results ────────────────────────────────────────────────────────────

let nextResultId = 0;
const queryResults = []; // [{ id, g1, g2, gm, fLabel }]

function freqLabel(freq) {
  return freq >= 1000
    ? `${(freq / 1000).toPrecision(4).replace(/\.?0+$/, '')} kHz`
    : `${freq} Hz`;
}

function addQueryResult(g1, g2, fLabel) {
  const gm = (g1 + g2) / 2;
  queryResults.push({ id: nextResultId++, g1, g2, gm, fLabel });
  renderQueryResults();
  // 末尾の新しい結果行にスクロール
  const list = document.getElementById('query-results-list');
  list.lastElementChild?.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}

function removeQueryResult(id) {
  const idx = queryResults.findIndex(r => r.id === id);
  if (idx !== -1) queryResults.splice(idx, 1);
  renderQueryResults();
}

function renderQueryResults() {
  const dp  = decimalPlaces;
  const fmt = v => (v >= 0 ? '+' : '') + v.toFixed(dp) + ' dB';
  const list    = document.getElementById('query-results-list');
  const section = document.getElementById('query-results-section');

  list.innerHTML = '';
  section.hidden = queryResults.length === 0;

  for (const res of queryResults) {
    const el = document.createElement('div');
    el.className = 'result-row';
    el.dataset.id = res.id;
    el.innerHTML = `
      <div class="result-cards">
        <div class="result-item result-item--eq1">
          <span class="result-label">EQ 1</span>
          <span class="result-value">${fmt(res.g1)}</span>
        </div>
        <div class="result-item result-item--eq2">
          <span class="result-label">EQ 2</span>
          <span class="result-value">${fmt(res.g2)}</span>
        </div>
        <div class="result-item result-item--mixed">
          <span class="result-label">ミックス @ ${res.fLabel}</span>
          <span class="result-value">${fmt(res.gm)}</span>
        </div>
      </div>
      <button class="btn-delete-row" aria-label="削除">
        <svg width="11" height="11" viewBox="0 0 12 12" fill="none">
          <path d="M1 1l10 10M11 1L1 11" stroke="currentColor" stroke-width="2" stroke-linecap="round"/>
        </svg>
      </button>
    `;
    list.appendChild(el);
  }
}

function doManualQuery() {
  const input = document.getElementById('query-freq-input');
  const freq  = parseFloat(input.value);
  if (!isFinite(freq) || freq <= 0) {
    input.classList.add('input-error');
    setTimeout(() => input.classList.remove('input-error'), 600);
    return;
  }
  const g1 = gainAtFreq(state.eq1, freq);
  const g2 = gainAtFreq(state.eq2, freq);
  addQueryResult(g1, g2, freqLabel(freq));
  input.value = '';
  input.focus();
}

function addQueryFromChart(rawFreq) {
  const rounded = Math.round(Math.max(FREQ_MIN, Math.min(FREQ_MAX, rawFreq)));
  const g1 = gainAtFreq(state.eq1, rounded);
  const g2 = gainAtFreq(state.eq2, rounded);
  addQueryResult(g1, g2, freqLabel(rounded));
}

// ─── File Format: FEQ ────────────────────────────────────────────────────────

/**
 * .feq テキストをパースして {freq, gain}[] を返す
 * フォーマット: 1行1ゲイン値、18バンド固定
 */
function parseFEQ(text) {
  const lines = text.trim().split(/\r?\n/);
  const gains = [];
  for (const line of lines) {
    const trimmed = line.trim();
    if (trimmed === '') continue;
    const val = parseFloat(trimmed);
    if (!isFinite(val)) throw new Error(`FEQ: 無効な値 "${trimmed}"`);
    gains.push(val);
  }
  if (gains.length !== 18) {
    throw new Error(`FEQ: 18バンド必要ですが ${gains.length} バンドが見つかりました`);
  }
  return FEQ_FREQS.map((freq, i) => ({ freq, gain: gains[i] }));
}

/**
 * ミックスカーブを .feq 形式テキストで生成
 * FEQ_FREQS の各周波数でミックスゲインを計算し、1行1値で出力
 */
function generateFEQ() {
  const lines = FEQ_FREQS.map(freq => {
    const g1 = gainAtFreq(state.eq1, freq);
    const g2 = gainAtFreq(state.eq2, freq);
    const gm = (g1 + g2) / 2;
    return gm.toFixed(1);
  });
  return lines.join('\n') + '\n';
}

// ─── File Format: XGEQ ───────────────────────────────────────────────────────

/**
 * .xgeq バイナリをパースして {freq, gain}[] を返す
 * 複数の戦略を試みて最良の結果を返す
 */
function parseXGEQ(buffer) {
  const bytes = new Uint8Array(buffer);

  // "v:" マーカーを探す (0x76, 0x3A)
  let dataStart = -1;
  for (let i = 0; i < bytes.length - 1; i++) {
    if (bytes[i] === 0x76 && bytes[i + 1] === 0x3A) {
      dataStart = i + 2;
      break;
    }
  }
  if (dataStart === -1) throw new Error('XGEQ: "v:" マーカーが見つかりません');

  // デバッグ: コンソールにhexダンプ
  const hexDump = Array.from(bytes.slice(dataStart, Math.min(dataStart + 96, bytes.length)))
    .map(b => b.toString(16).padStart(2, '0')).join(' ');
  console.log('[XGEQ] total bytes:', buffer.byteLength, '/ after v::', buffer.byteLength - dataStart);
  console.log('[XGEQ] hex after v::', hexDump);

  const view = new DataView(buffer);
  let bestBands = [];

  const isValidFreq = f => isFinite(f) && f >= 10 && f <= 30000;
  const isValidGain = g => isFinite(g) && Math.abs(g) <= 120;
  const makeBand    = (f, g) => ({ freq: Math.round(f * 10) / 10, gain: Math.round(g * 100) / 100 });
  const tryUpdate   = arr => { if (arr.length > bestBands.length) bestBands = arr; };

  // ── 戦略 A: float32 ペア (freq, gain) ── オフセット 0-20、LE/BE ──
  for (let extra = 0; extra <= 20; extra++) {
    const off = dataStart + extra;
    if (off + 8 > buffer.byteLength) break;
    for (const le of [true, false]) {
      const n = Math.floor((buffer.byteLength - off) / 8);
      const bands = [];
      for (let i = 0; i < n; i++) {
        const freq = view.getFloat32(off + i * 8,     le);
        const gain = view.getFloat32(off + i * 8 + 4, le);
        if (isValidFreq(freq) && isValidGain(gain)) bands.push(makeBand(freq, gain));
      }
      if (bands.length > 0) {
        console.log(`[XGEQ] A float32-pairs extra=${extra} le=${le}:`, bands.length, 'bands', bands[0]);
        tryUpdate(bands);
      }
    }
  }

  // ── 戦略 B: uint32カウント + float32 ペア ── オフセット 0-12、LE/BE ──
  for (let extra = 0; extra <= 12; extra++) {
    const off = dataStart + extra;
    if (off + 4 > buffer.byteLength) break;
    for (const le of [true, false]) {
      const count = view.getUint32(off, le);
      if (count > 0 && count <= 512 && off + 4 + count * 8 <= buffer.byteLength) {
        const bands = [];
        for (let i = 0; i < count; i++) {
          const freq = view.getFloat32(off + 4 + i * 8,     le);
          const gain = view.getFloat32(off + 4 + i * 8 + 4, le);
          if (isValidFreq(freq) && isValidGain(gain)) bands.push(makeBand(freq, gain));
        }
        if (bands.length > 0) {
          console.log(`[XGEQ] B count+float32 extra=${extra} le=${le} count=${count}:`, bands.length);
          tryUpdate(bands);
        }
      }
    }
  }

  // ── 戦略 C: float64 ペア ── オフセット 0-8、LE ──
  for (let extra = 0; extra <= 8; extra++) {
    const off = dataStart + extra;
    if (off + 16 > buffer.byteLength) break;
    const n = Math.floor((buffer.byteLength - off) / 16);
    const bands = [];
    for (let i = 0; i < n; i++) {
      const freq = view.getFloat64(off + i * 16,     true);
      const gain = view.getFloat64(off + i * 16 + 8, true);
      if (isValidFreq(freq) && isValidGain(gain)) bands.push(makeBand(freq, gain));
    }
    if (bands.length > 0) {
      console.log(`[XGEQ] C float64-pairs extra=${extra}:`, bands.length, bands[0]);
      tryUpdate(bands);
    }
  }

  // ── 戦略 D: uint32カウント + float64 ペア ─────────────────────────────
  for (let extra = 0; extra <= 8; extra++) {
    const off = dataStart + extra;
    if (off + 4 > buffer.byteLength) break;
    for (const le of [true, false]) {
      const count = view.getUint32(off, le);
      if (count > 0 && count <= 256 && off + 4 + count * 16 <= buffer.byteLength) {
        const bands = [];
        for (let i = 0; i < count; i++) {
          const freq = view.getFloat64(off + 4 + i * 16,     le);
          const gain = view.getFloat64(off + 4 + i * 16 + 8, le);
          if (isValidFreq(freq) && isValidGain(gain)) bands.push(makeBand(freq, gain));
        }
        if (bands.length > 0) {
          console.log(`[XGEQ] D count+float64 extra=${extra} le=${le} count=${count}:`, bands.length);
          tryUpdate(bands);
        }
      }
    }
  }

  // ── 戦略 E: ゲインのみ float32 (固定周波数 FEQ_FREQS 18バンド) ─────────
  for (let extra = 0; extra <= 12; extra++) {
    const off = dataStart + extra;
    if (off + 18 * 4 > buffer.byteLength) break;
    const bands = [];
    for (let i = 0; i < 18; i++) {
      const gain = view.getFloat32(off + i * 4, true);
      if (isValidGain(gain)) bands.push({ freq: FEQ_FREQS[i], gain: Math.round(gain * 100) / 100 });
    }
    if (bands.length === 18) {
      console.log(`[XGEQ] E gains-only-18 extra=${extra}:`, bands.slice(0, 3));
      tryUpdate(bands);
    }
  }

  if (bestBands.length > 0) return bestBands;
  throw new Error('XGEQ: 有効なバンドデータが見つかりません（F12キー→コンソールタブで [XGEQ] を確認してください）');
}

/**
 * ミックスカーブを .xgeq バイナリで生成
 * フォーマット: "foo_dsp_xgeq\n1\nv:" + uint32LE(count) + N*(float32LE freq, float32LE gain)
 */
function generateXGEQ() {
  const freqs = XGEQ_EXPORT_FREQS;
  const n = freqs.length;

  const headerStr = 'foo_dsp_xgeq\n1\nv:';
  const headerBytes = new TextEncoder().encode(headerStr);

  // uint32 count + n * (float32 freq + float32 gain)
  const binBuffer = new ArrayBuffer(4 + n * 8);
  const view = new DataView(binBuffer);
  view.setUint32(0, n, true); // LE
  for (let i = 0; i < n; i++) {
    const freq = freqs[i];
    const g1   = gainAtFreq(state.eq1, freq);
    const g2   = gainAtFreq(state.eq2, freq);
    const gain = (g1 + g2) / 2;
    view.setFloat32(4 + i * 8,     freq, true);
    view.setFloat32(4 + i * 8 + 4, gain, true);
  }

  const combined = new Uint8Array(headerBytes.length + binBuffer.byteLength);
  combined.set(headerBytes, 0);
  combined.set(new Uint8Array(binBuffer), headerBytes.length);
  return combined;
}

// ─── File Download Helper ─────────────────────────────────────────────────────

function downloadFile(filename, data) {
  const blob = data instanceof Uint8Array
    ? new Blob([data], { type: 'application/octet-stream' })
    : new Blob([data], { type: 'text/plain; charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

// ─── Modal ────────────────────────────────────────────────────────────────────

const modalState = {
  mode:   null, // 'import' | 'export'
  target: null, // 'eq1' | 'eq2' (importのみ)
};

function openModal(mode, target) {
  modalState.mode   = mode;
  modalState.target = target || null;

  const overlay    = document.getElementById('modal-overlay');
  const title      = document.getElementById('modal-title');
  const fileField  = document.getElementById('modal-file-field');
  const fileInput  = document.getElementById('modal-file-input');
  const errorEl    = document.getElementById('modal-error');
  const confirmBtn = document.getElementById('modal-confirm');

  // タイトル
  if (mode === 'import') {
    title.textContent = `インポート（${target === 'eq1' ? 'EQ 1' : 'EQ 2'}）`;
    fileField.hidden  = false;
    confirmBtn.textContent = 'インポート';
  } else {
    title.textContent = 'エクスポート（ミックスカーブ）';
    fileField.hidden  = true;
    confirmBtn.textContent = 'ダウンロード';
  }

  // リセット
  fileInput.value = '';
  errorEl.hidden  = true;
  errorEl.textContent = '';

  // format ラジオを feq に戻す
  const radios = overlay.querySelectorAll('input[name="modal-format"]');
  radios.forEach(r => { r.checked = r.value === 'feq'; });

  overlay.hidden = false;
  // アクセシビリティ: フォーカスをモーダルに移動
  document.getElementById('modal-close').focus();
}

function closeModal() {
  document.getElementById('modal-overlay').hidden = true;
  modalState.mode   = null;
  modalState.target = null;
}

function getSelectedFormat() {
  const checked = document.querySelector('input[name="modal-format"]:checked');
  return checked ? checked.value : 'feq';
}

function showModalError(msg) {
  const errorEl = document.getElementById('modal-error');
  errorEl.textContent = msg;
  errorEl.hidden = false;
}

function handleModalConfirm() {
  const format = getSelectedFormat();

  if (modalState.mode === 'export') {
    // ─ エクスポート ─
    try {
      if (format === 'feq') {
        const text = generateFEQ();
        downloadFile('mixed_eq.feq', text);
      } else {
        const bytes = generateXGEQ();
        downloadFile('mixed_eq.xgeq', bytes);
      }
      closeModal();
    } catch (e) {
      showModalError(`エクスポートエラー: ${e.message}`);
    }

  } else if (modalState.mode === 'import') {
    // ─ インポート ─
    const fileInput = document.getElementById('modal-file-input');
    const file = fileInput.files[0];
    if (!file) {
      showModalError('ファイルを選択してください');
      return;
    }

    const eqKey = modalState.target;

    if (format === 'feq') {
      const reader = new FileReader();
      reader.onload = e => {
        try {
          const bands = parseFEQ(e.target.result);
          state[eqKey] = bands;
          renderBands(eqKey);
          updateChart();
          closeModal();
        } catch (err) {
          showModalError(err.message);
        }
      };
      reader.onerror = () => showModalError('ファイルの読み込みに失敗しました');
      reader.readAsText(file);

    } else {
      const reader = new FileReader();
      reader.onload = e => {
        try {
          const bands = parseXGEQ(e.target.result);
          state[eqKey] = bands;
          renderBands(eqKey);
          updateChart();
          closeModal();
        } catch (err) {
          showModalError(err.message);
        }
      };
      reader.onerror = () => showModalError('ファイルの読み込みに失敗しました');
      reader.readAsArrayBuffer(file);
    }
  }
}

// ─── Boot ─────────────────────────────────────────────────────────────────────

document.addEventListener('DOMContentLoaded', () => {
  renderBands('eq1');
  renderBands('eq2');

  // イベントリスナーを先に設定し、chart初期化エラーの影響を受けないようにする

  // Add band
  document.getElementById('add-eq1').addEventListener('click', () => addBand('eq1'));
  document.getElementById('add-eq2').addEventListener('click', () => addBand('eq2'));

  // Clear
  document.getElementById('clear-eq1').addEventListener('click', () => clearEQ('eq1'));
  document.getElementById('clear-eq2').addEventListener('click', () => clearEQ('eq2'));

  // Import
  document.getElementById('import-eq1').addEventListener('click', () => openModal('import', 'eq1'));
  document.getElementById('import-eq2').addEventListener('click', () => openModal('import', 'eq2'));

  // Export
  document.getElementById('export-btn').addEventListener('click', () => openModal('export'));

  // Input / remove — event delegation per band-list container
  for (const eqKey of ['eq1', 'eq2']) {
    const container = document.getElementById(`bands-${eqKey}`);

    container.addEventListener('input', e => {
      const input = e.target.closest('.band-input');
      if (!input) return;
      updateBand(input.dataset.eq, Number(input.dataset.idx), input.dataset.field, input.value);
    });

    container.addEventListener('click', e => {
      const btn = e.target.closest('.btn-remove');
      if (!btn) return;
      removeBand(btn.dataset.eq, Number(btn.dataset.idx));
    });
  }

  // Query: 取得ボタン・Enterキー
  document.getElementById('query-submit-btn').addEventListener('click', doManualQuery);
  document.getElementById('query-freq-input').addEventListener('keydown', e => {
    if (e.key === 'Enter') doManualQuery();
  });

  // Query: 結果の削除（イベント委譲）
  document.getElementById('query-results-list').addEventListener('click', e => {
    const btn = e.target.closest('.btn-delete-row');
    if (!btn) return;
    const row = btn.closest('.result-row');
    if (row) removeQueryResult(Number(row.dataset.id));
  });

  // 小数点桁数
  document.getElementById('decimal-places').addEventListener('change', e => {
    decimalPlaces = Number(e.target.value);
    renderQueryResults();
  });

  // Modal: 閉じる / キャンセル / 確定 / オーバーレイクリック
  document.getElementById('modal-close').addEventListener('click', closeModal);
  document.getElementById('modal-cancel').addEventListener('click', closeModal);
  document.getElementById('modal-confirm').addEventListener('click', handleModalConfirm);
  document.getElementById('modal-overlay').addEventListener('click', e => {
    if (e.target === e.currentTarget) closeModal();
  });
  // Escape キーで閉じる
  document.addEventListener('keydown', e => {
    if (e.key === 'Escape' && !document.getElementById('modal-overlay').hidden) {
      closeModal();
    }
  });

  // format ラジオ切替でファイル入力の accept を更新
  document.querySelectorAll('input[name="modal-format"]').forEach(radio => {
    radio.addEventListener('change', () => {
      const fmt = getSelectedFormat();
      const fileInput = document.getElementById('modal-file-input');
      fileInput.accept = fmt === 'feq' ? '.feq' : '.xgeq';
      // エラー表示をリセット
      const errorEl = document.getElementById('modal-error');
      errorEl.hidden = true;
    });
  });

  // Chartの初期化はイベントリスナー設定後に行う
  try {
    initChart();
    updateChart();
  } catch (e) {
    console.error('Chart initialization failed:', e);
  }
});
