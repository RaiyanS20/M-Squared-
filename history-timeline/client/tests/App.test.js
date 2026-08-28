import { describe, it, before, beforeEach } from 'node:test';
import assert from 'node:assert/strict';
import { setupDom, APP_HTML, appRoots, flush } from './domEnvironment.js';
import { COUNTRIES_1969, EVENTS_FRA_1969, TIMELINE } from './stubApi.js';

before(() => setupDom(APP_HTML));

const { App } = await import('../js/App.js');
const { TimelineApi } = await import('../js/api/TimelineApi.js');
const { stubFetch } = await import('./stubApi.js');

/**
 * THE WHOLE APP, driven the way a learner uses it, with only the network
 * stubbed. This one file proves all four MVP features work TOGETHER — which is
 * the assertion that actually matters to the product.
 *
 * These are few and slow(ish) by design: when a whole-app test fails you know
 * *something* broke; the unit tests tell you *what*. You want many of those and
 * a handful of these.
 */
describe('App', () => {
  let roots;
  let stub;

  const happyPath = () =>
    stubFetch({
      // Order matters: stubFetch matches by substring, so the longer, more
      // specific path must be listed before the prefix it contains.
      '/api/years/1969/countries/FRA/events': EVENTS_FRA_1969,
      '/api/years/1969/countries': COUNTRIES_1969,
      '/api/timeline': TIMELINE,
    });

  const text = (key) => roots[key].textContent;
  const button = (label) =>
    [...roots.countries.querySelectorAll('button')].find((b) => b.textContent.includes(label));
  const year = (y) => roots.timeline.querySelector(`[data-year="${y}"]`);
  const checkedYear = () =>
    roots.timeline.querySelector('[role="radio"][aria-checked="true"]')?.dataset.year;

  beforeEach(() => {
    document.body.innerHTML = APP_HTML;
    roots = appRoots();
  });

  async function startApp() {
    const app = new App({ api: new TimelineApi(), roots });
    await app.start();
    await flush();
    return app;
  }

  it('walks a learner from year to country to events', async () => {
    stub = happyPath();
    await startApp();

    // Feature 1 + the sensible default: the most recent year that actually has
    // events (1969), not the empty 1970.
    assert.equal(checkedYear(), '1969');

    // Feature 3: countries for that year.
    assert.ok(button('France'), 'expected France in the country list');
    assert.ok(button('Japan'));

    // Feature 4: pick one and read what happened.
    button('France').click();
    await flush();
    assert.match(text('events'), /Concorde first flight/);
    assert.match(text('events'), /De Gaulle resigns/);

    stub.restore();
  });

  it('prompts for a country before any is chosen', async () => {
    stub = happyPath();
    await startApp();
    assert.match(text('events'), /choose a country to see its events/i);
    stub.restore();
  });

  /**
   * THE STALE-DATA TEST.
   *
   * Changing the year clears the country. If the events slice were not reset to
   * `idle` at the same time, the PREVIOUS country's events would stay on screen
   * under the new year's heading — plausible-looking and completely wrong, and
   * invisible unless you happen to look closely. This is what tests are for.
   */
  it('clears the chosen country and its events when the year changes', async () => {
    stub = stubFetch({
      '/api/years/1969/countries/FRA/events': EVENTS_FRA_1969,
      '/api/years/1969/countries': COUNTRIES_1969,
      '/api/years/1967/countries': { year: 1967, countries: [] },
      '/api/timeline': TIMELINE,
    });
    const app = await startApp();

    button('France').click();
    await flush();
    assert.match(text('events'), /Concorde first flight/);

    // Feature 2: move to a different year.
    year(1967).click();
    await flush();

    assert.doesNotMatch(text('events'), /Concorde first flight/);
    assert.match(text('events'), /choose a country to see its events/i);
    assert.equal(app.state.selectedCountryCode, null);

    stub.restore();
  });

  it('explains a year with nothing recorded instead of showing a blank panel', async () => {
    stub = stubFetch({
      '/api/years/1967/countries': { year: 1967, countries: [] },
      '/api/years/1969/countries': COUNTRIES_1969,
      '/api/timeline': TIMELINE,
    });
    await startApp();

    year(1967).click();
    await flush();
    assert.match(text('countries'), /nothing is recorded for 1967/i);

    stub.restore();
  });

  it('surfaces a server error instead of failing silently', async () => {
    stub = stubFetch(
      { '/api/timeline': { error: { code: 'INTERNAL_ERROR', message: 'The archive is offline.' } } },
      { status: { '/api/timeline': 500 } },
    );
    await startApp();

    assert.match(text('timeline'), /the archive is offline/i);
    assert.ok(roots.timeline.querySelector('[role="alert"]'), 'errors must be announced');

    stub.restore();
  });

  it('ignores a repeated click on the year already selected', async () => {
    stub = happyPath();
    await startApp();
    const before = stub.calls.length;

    year(1969).click();
    await flush();

    assert.equal(stub.calls.length, before, 'should not re-request the same year');
    stub.restore();
  });

  it('keeps the state minimal — only the two user choices are stored', async () => {
    stub = happyPath();
    const app = await startApp();
    button('France').click();
    await flush();

    // Everything else on screen is DERIVED from these two values.
    assert.equal(app.state.selectedYear, 1969);
    assert.equal(app.state.selectedCountryCode, 'FRA');
    assert.equal(app.state.events.status, 'success');

    stub.restore();
  });
});
