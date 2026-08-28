import { JSDOM } from 'jsdom';

/**
 * Gives Node a DOM.
 *
 * Node has no `document` — it is a server runtime. jsdom implements the DOM in
 * pure JavaScript so browser code can be tested without launching a browser.
 * That is the entire trick behind every front-end test runner you will meet;
 * doing it explicitly here means there is no configuration file doing something
 * invisible on your behalf.
 *
 * The trade-off is worth knowing: jsdom is fast and scriptable but it does not
 * lay anything out, so it has no notion of what is visible, scrolled or
 * overlapping. Tests that depend on real layout need a real browser
 * (Playwright); everything else is quicker and steadier here.
 */
export function setupDom(bodyHtml = '') {
  const dom = new JSDOM(`<!doctype html><html><body>${bodyHtml}</body></html>`, {
    url: 'http://localhost/',
    pretendToBeVisual: true,
  });

  // Copy the handful of globals browser code expects. Doing this by hand is
  // exactly what a test framework's "jsdom environment" does for you.
  //
  // `defineProperty` rather than plain assignment, because some of these names
  // (`navigator`) are read-only built-ins in Node, and assigning to a read-only
  // global throws in a module — modules are always in strict mode.
  const globals = ['window', 'document', 'Node', 'HTMLElement', 'Event', 'KeyboardEvent', 'MouseEvent', 'CustomEvent', 'DOMException', 'navigator'];
  for (const name of globals) {
    Object.defineProperty(globalThis, name, {
      value: dom.window[name],
      writable: true,
      configurable: true,
    });
  }

  // jsdom does not implement scrollIntoView (it does not lay anything out), and
  // TimelineView calls it. A no-op stub is enough.
  dom.window.Element.prototype.scrollIntoView = () => {};

  return dom;
}

/** The mount points that index.html provides. */
export const APP_HTML = `
  <div id="timeline"></div>
  <div id="countries"></div>
  <span id="countries-context"></span>
  <div id="events"></div>
`;

export function appRoots() {
  return {
    timeline: document.getElementById('timeline'),
    countries: document.getElementById('countries'),
    countriesContext: document.getElementById('countries-context'),
    events: document.getElementById('events'),
  };
}

/**
 * Lets pending promises settle.
 *
 * The app loads data in `async` methods, so after a click the DOM is not
 * updated until those promises resolve. `setTimeout(0)` yields to the event
 * loop, which is enough for promise chains that are not waiting on real I/O.
 */
export function flush(times = 3) {
  return new Promise((resolve) => {
    let remaining = times;
    const tick = () => (remaining-- > 0 ? setTimeout(tick, 0) : resolve());
    tick();
  });
}
