import { describe, it, before, beforeEach } from 'node:test';
import assert from 'node:assert/strict';
import { setupDom } from './domEnvironment.js';
import { COUNTRIES_1969, EVENTS_FRA_1969, FRANCE, JAPAN, TIMELINE, anEvent } from './stubApi.js';

before(() => setupDom());

const { TimelineView } = await import('../js/components/TimelineView.js');
const { CountryListView, groupByRegion } = await import('../js/components/CountryListView.js');
const { EventListView, formatMonthDay } = await import('../js/components/EventListView.js');
const { Async } = await import('../js/Store.js');

/**
 * View tests written the way a USER experiences the view: find things by their
 * accessible role or visible text, click them, assert on what appears.
 *
 * Almost nothing here asserts on a CSS class, and where it does (the selected
 * year) the class is the only way to observe a visual state. The heuristic:
 * would this test still pass if I rewrote the internals but kept the behaviour?
 * If no, you are testing the wrong thing.
 */

function mount() {
  const root = document.createElement('div');
  document.body.replaceChildren(root);
  return root;
}

const state = (over = {}) => ({
  selectedYear: null,
  selectedCountryCode: null,
  timeline: Async.idle(),
  countries: Async.idle(),
  events: Async.idle(),
  ...over,
});

describe('TimelineView', () => {
  let root;
  let selected;
  let view;

  beforeEach(() => {
    root = mount();
    selected = [];
    view = new TimelineView(root, { onSelectYear: (y) => selected.push(y) });
  });

  const ready = (year = 1969) =>
    state({ timeline: Async.success(TIMELINE), selectedYear: year });

  it('shows a loading message before the data arrives', () => {
    view.update(state({ timeline: Async.loading() }));
    assert.match(root.textContent, /loading the timeline/i);
  });

  it('shows the error message when the request fails', () => {
    view.update(state({ timeline: Async.error({ message: 'The archive is offline.' }) }));
    assert.match(root.textContent, /archive is offline/i);
    assert.equal(root.querySelector('[role="alert"]') !== null, true);
  });

  it('renders one radio per year inside a radiogroup', () => {
    view.update(ready());
    assert.equal(root.querySelector('[role="radiogroup"]') !== null, true);
    assert.equal(root.querySelectorAll('[role="radio"]').length, 4);
  });

  it('labels each year with its event count for screen readers', () => {
    view.update(ready());
    const y1969 = root.querySelector('[data-year="1969"]');
    const y1968 = root.querySelector('[data-year="1968"]');
    assert.equal(y1969.getAttribute('aria-label'), '1969, 2 events');
    assert.equal(y1968.getAttribute('aria-label'), '1968, 1 event'); // singular
  });

  it('marks exactly one year as checked', () => {
    view.update(ready(1969));
    const checked = [...root.querySelectorAll('[role="radio"]')].filter(
      (r) => r.getAttribute('aria-checked') === 'true',
    );
    assert.equal(checked.length, 1);
    assert.equal(checked[0].dataset.year, '1969');
  });

  it('keeps exactly one year in the tab order (roving tabindex)', () => {
    view.update(ready(1969));
    const tabbable = [...root.querySelectorAll('[role="radio"]')].filter(
      (r) => r.getAttribute('tabindex') === '0',
    );
    assert.equal(tabbable.length, 1);
  });

  it('reports the year that was clicked', () => {
    view.update(ready());
    root.querySelector('[data-year="1967"]').click();
    assert.deepEqual(selected, [1967]);
  });

  /**
   * This is the test that protects the optimisation described at the top of
   * TimelineView: changing the SELECTION must not rebuild the buttons, or the
   * focused element would be destroyed mid-keypress.
   */
  it('reuses the same button elements when only the selection changes', () => {
    view.update(ready(1969));
    const before = root.querySelector('[data-year="1968"]');
    view.update(ready(1968));
    const after = root.querySelector('[data-year="1968"]');
    assert.equal(before, after, 'the button should be the very same DOM node');
    assert.equal(after.getAttribute('aria-checked'), 'true');
  });

  it('rebuilds when the timeline data itself changes', () => {
    view.update(ready(1969));
    const before = root.querySelector('[data-year="1969"]');
    const newTimeline = { ...TIMELINE, ticks: TIMELINE.ticks.map((t) => ({ ...t })) };
    view.update(state({ timeline: Async.success(newTimeline), selectedYear: 1969 }));
    assert.notEqual(before, root.querySelector('[data-year="1969"]'));
  });

  describe('keyboard navigation', () => {
    const press = (key) => {
      const event = new window.KeyboardEvent('keydown', { key, bubbles: true, cancelable: true });
      root.querySelector('[role="radiogroup"]').dispatchEvent(event);
    };

    beforeEach(() => view.update(ready(1969)));

    it('moves one year with the arrow keys', () => {
      press('ArrowRight');
      assert.deepEqual(selected, [1970]);
      press('ArrowLeft');
      assert.equal(selected.at(-1), 1968);
    });

    it('jumps a decade with Page Up and Page Down, clamped to the range', () => {
      press('PageDown');
      assert.equal(selected.at(-1), 1970); // clamped at the last year
      press('PageUp');
      assert.equal(selected.at(-1), 1967); // clamped at the first
    });

    it('jumps to the ends with Home and End', () => {
      press('Home');
      assert.equal(selected.at(-1), 1967);
      press('End');
      assert.equal(selected.at(-1), 1970);
    });

    it('clamps instead of running off the end', () => {
      view.update(ready(1970));
      press('ArrowRight');
      assert.equal(selected.at(-1), 1970);
    });

    it('ignores keys it does not handle', () => {
      press('a');
      assert.deepEqual(selected, []);
    });
  });
});

describe('groupByRegion', () => {
  it('groups countries under their region, sorted', () => {
    const grouped = groupByRegion([JAPAN, FRANCE]);
    assert.deepEqual(grouped.map(([region]) => region), ['Asia', 'Europe']);
    assert.deepEqual(grouped[1][1].map((c) => c.name), ['France']);
  });

  it('sorts countries within a region by name', () => {
    const spain = { id: 9, code: 'ESP', name: 'Spain', region: 'Europe' };
    const grouped = groupByRegion([spain, FRANCE]);
    assert.deepEqual(grouped[0][1].map((c) => c.name), ['France', 'Spain']);
  });

  it('returns an empty array for no countries', () => {
    assert.deepEqual(groupByRegion([]), []);
  });
});

describe('CountryListView', () => {
  let root;
  let context;
  let chosen;
  let view;

  beforeEach(() => {
    root = mount();
    context = document.createElement('span');
    chosen = [];
    view = new CountryListView(root, {
      contextRoot: context,
      onSelectCountry: (code) => chosen.push(code),
    });
  });

  it('prompts for a year before one is chosen', () => {
    view.update(state());
    assert.match(root.textContent, /select a year to begin/i);
    assert.equal(context.textContent, '');
  });

  it('shows the chosen year beside the heading', () => {
    view.update(state({ selectedYear: 1969, countries: Async.success(COUNTRIES_1969) }));
    assert.equal(context.textContent, ' in 1969');
  });

  it('groups countries by region', () => {
    view.update(state({ selectedYear: 1969, countries: Async.success(COUNTRIES_1969) }));
    const headings = [...root.querySelectorAll('h3')].map((h) => h.textContent);
    assert.deepEqual(headings, ['Asia', 'Europe']);
  });

  it('reports the chosen country code', () => {
    view.update(state({ selectedYear: 1969, countries: Async.success(COUNTRIES_1969) }));
    [...root.querySelectorAll('button')].find((b) => b.textContent.includes('Japan')).click();
    assert.deepEqual(chosen, ['JPN']);
  });

  it('marks the selected country for assistive technology', () => {
    view.update(
      state({ selectedYear: 1969, selectedCountryCode: 'FRA', countries: Async.success(COUNTRIES_1969) }),
    );
    const france = [...root.querySelectorAll('button')].find((b) => b.textContent.includes('France'));
    assert.equal(france.getAttribute('aria-current'), 'true');
  });

  it('explains a year with nothing recorded rather than showing a blank panel', () => {
    view.update(state({ selectedYear: 1902, countries: Async.success({ year: 1902, countries: [] }) }));
    assert.match(root.textContent, /nothing is recorded for 1902/i);
  });
});

describe('formatMonthDay', () => {
  it('formats a valid month and day', () => {
    assert.equal(formatMonthDay('03-02'), '2 March');
    assert.equal(formatMonthDay('12-25'), '25 December');
    assert.equal(formatMonthDay('01-01'), '1 January');
  });

  it('returns null for anything it cannot format', () => {
    for (const bad of [null, '', '13-01', '00-05', 'nonsense']) {
      assert.equal(formatMonthDay(bad), null, `should reject ${JSON.stringify(bad)}`);
    }
  });
});

describe('EventListView', () => {
  let root;
  let view;

  beforeEach(() => {
    root = mount();
    view = new EventListView(root);
  });

  it('prompts for a country before one is chosen', () => {
    view.update(state());
    assert.match(root.textContent, /choose a country to see its events/i);
  });

  it('renders each event with its date, summary and category', () => {
    view.update(state({ events: Async.success(EVENTS_FRA_1969) }));
    assert.match(root.textContent, /Concorde first flight/);
    assert.match(root.textContent, /2 March 1969/);
    assert.match(root.textContent, /science/);
    assert.equal(root.querySelectorAll('li').length, 2);
  });

  it('shows just the year when the exact date is unknown', () => {
    view.update(
      state({ events: Async.success({ ...EVENTS_FRA_1969, events: [anEvent({ monthDay: null })] }) }),
    );
    assert.equal(root.querySelector('time').textContent, '1969');
  });

  it('exposes a machine-readable date', () => {
    view.update(state({ events: Async.success(EVENTS_FRA_1969) }));
    assert.equal(root.querySelector('time').getAttribute('datetime'), '1969-03-02');
  });

  // The editorial promise of the project, verified in the UI.
  it('gives every event a source link that opens safely in a new tab', () => {
    view.update(state({ events: Async.success(EVENTS_FRA_1969) }));
    const links = [...root.querySelectorAll('a')];
    assert.equal(links.length, EVENTS_FRA_1969.events.length);
    for (const link of links) {
      assert.match(link.getAttribute('href'), /^https?:\/\//);
      assert.equal(link.getAttribute('rel'), 'noopener noreferrer');
    }
  });

  it('explains an empty result rather than showing a blank panel', () => {
    view.update(state({ events: Async.success({ ...EVENTS_FRA_1969, events: [] }) }));
    assert.match(root.textContent, /nothing is recorded for France in 1969/i);
  });

  // Untrusted content reaching the DOM is the risk this whole helper guards
  // against, so it is asserted at the view level too.
  it('renders a malicious title as text, never as markup', () => {
    const attack = '<img src=x onerror="globalThis.__pwnedView = true">';
    view.update(
      state({ events: Async.success({ ...EVENTS_FRA_1969, events: [anEvent({ title: attack })] }) }),
    );
    assert.equal(root.querySelector('img'), null);
    assert.equal(globalThis.__pwnedView, undefined);
    assert.match(root.textContent, /onerror/);
  });
});
