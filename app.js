'use strict';

// ─── Constants ────────────────────────────────────────────────────────────────

const FREQ_MIN   = 20;
const FREQ_MAX   = 20000;
const N_POINTS   = 400;

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

// ─── Query Rows ───────────────────────────────────────────────────────────────

let nextQueryId = 1;
// [{ id, freq, result: null | { g1, g2, gm, fLabel } }]
// id=0 は index.html に静的に存在する初期行と対応
const queryRows = [{ id: 0, freq: '', result: null }];

function freqLabel(freq) {
  return freq >= 1000
    ? `${(freq / 1000).toPrecision(4).replace(/\.?0+$/, '')} kHz`
    : `${freq} Hz`;
}

function addQueryFromChart(rawFreq) {
  const freq = Math.max(FREQ_MIN, Math.min(FREQ_MAX, rawFreq));
  // 有効数字3桁に丸める（例: 1234.5 → 1230, 45.6 → 45.6）
  const magnitude = Math.pow(10, Math.floor(Math.log10(freq)) - 2);
  const rounded = Math.round(freq / magnitude) * magnitude;

  const id = nextQueryId++;
  const g1 = gainAtFreq(state.eq1, rounded);
  const g2 = gainAtFreq(state.eq2, rounded);
  queryRows.push({
    id,
    freq: String(rounded),
    result: { g1, g2, gm: (g1 + g2) / 2, fLabel: freqLabel(rounded) },
  });
  renderQueryRows();

  const el = document.querySelector(`.query-row[data-id="${id}"]`);
  if (el) el.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}

function addQueryRow() {
  const id = nextQueryId++;
  queryRows.push({ id, freq: '', result: null });
  renderQueryRows();
  // Focus the new input
  const input = document.querySelector(`.query-row[data-id="${id}"] .query-input`);
  if (input) input.focus();
}

function removeQueryRow(id) {
  const idx = queryRows.findIndex(r => r.id === id);
  if (idx === -1) return;
  queryRows.splice(idx, 1);
  if (queryRows.length === 0) addQueryRow();
  else renderQueryRows();
}

function executeQuery(id) {
  const row  = queryRows.find(r => r.id === id);
  if (!row) return;

  const freq = parseFloat(row.freq);
  if (!isFinite(freq) || freq <= 0) {
    row.result = { error: '有効な周波数を入力してください（例: 1000）' };
    renderQueryRows();
    return;
  }

  const g1 = gainAtFreq(state.eq1, freq);
  const g2 = gainAtFreq(state.eq2, freq);
  row.result = { g1, g2, gm: (g1 + g2) / 2, fLabel: freqLabel(freq) };
  renderQueryRows();

  // Add new empty row only if this is the last one
  const isLast = queryRows[queryRows.length - 1].id === id;
  if (isLast) addQueryRow();
}

function renderQueryRows() {
  const fmt = v => (v >= 0 ? '+' : '') + v.toFixed(2) + ' dB';
  const container = document.getElementById('query-rows');
  container.innerHTML = '';

  for (const row of queryRows) {
    const el = document.createElement('div');
    el.className = 'query-row';
    el.dataset.id = row.id;

    let resultHTML = '';
    if (row.result) {
      if (row.result.error) {
        resultHTML = `<div class="query-row-result"><span class="result-error">${row.result.error}</span></div>`;
      } else {
        const { g1, g2, gm, fLabel } = row.result;
        resultHTML = `
          <div class="query-row-result">
            <div class="result-grid">
              <div class="result-item result-item--eq1">
                <span class="result-label">EQ 1</span>
                <span class="result-value">${fmt(g1)}</span>
              </div>
              <div class="result-item result-item--eq2">
                <span class="result-label">EQ 2</span>
                <span class="result-value">${fmt(g2)}</span>
              </div>
              <div class="result-item result-item--mixed">
                <span class="result-label">ミックス @ ${fLabel}</span>
                <span class="result-value">${fmt(gm)}</span>
              </div>
            </div>
          </div>`;
      }
    }

    el.innerHTML = `
      <div class="query-row-form">
        <div class="query-input-wrap">
          <input
            type="number"
            class="query-input"
            placeholder="1000"
            min="1" max="96000" step="1"
            value="${row.freq}"
            aria-label="周波数 (Hz)"
          >
          <span class="query-unit">Hz</span>
        </div>
        <button class="btn-primary btn-get">取得</button>
        ${row.result ? `
        <button class="btn-delete-row" aria-label="この行を削除">
          <svg width="11" height="11" viewBox="0 0 12 12" fill="none">
            <path d="M1 1l10 10M11 1L1 11" stroke="currentColor" stroke-width="2" stroke-linecap="round"/>
          </svg>
          削除
        </button>` : ''}
      </div>
      ${resultHTML}
    `;

    container.appendChild(el);
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

  // Query rows — event delegation on container
  document.getElementById('query-rows').addEventListener('click', e => {
    const row = e.target.closest('.query-row');
    if (!row) return;
    const id = Number(row.dataset.id);

    if (e.target.closest('.btn-get')) {
      executeQuery(id);
    } else if (e.target.closest('.btn-delete-row')) {
      removeQueryRow(id);
    }
  });

  document.getElementById('query-rows').addEventListener('keydown', e => {
    if (e.key !== 'Enter') return;
    const input = e.target.closest('.query-input');
    if (!input) return;
    const row = input.closest('.query-row');
    if (!row) return;
    executeQuery(Number(row.dataset.id));
  });

  document.getElementById('query-rows').addEventListener('input', e => {
    const input = e.target.closest('.query-input');
    if (!input) return;
    const row = input.closest('.query-row');
    if (!row) return;
    const r = queryRows.find(r => r.id === Number(row.dataset.id));
    if (r) r.freq = input.value;
  });

  // Chartの初期化はイベントリスナー設定後に行う
  try {
    initChart();
    updateChart();
  } catch (e) {
    console.error('Chart initialization failed:', e);
  }
});
