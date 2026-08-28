import { App } from './App.js';
import { timelineApi } from './api/TimelineApi.js';

/**
 * The entry point — the browser equivalent of the server's composition root.
 *
 * Find the mount points, build the app with its real dependencies, start it.
 * Nothing else. Keeping this separate from `App.js` is what lets the tests
 * construct an `App` with a stubbed API and their own DOM.
 */
const roots = {
  timeline: document.getElementById('timeline'),
  countries: document.getElementById('countries'),
  countriesContext: document.getElementById('countries-context'),
  events: document.getElementById('events'),
};

for (const [name, node] of Object.entries(roots)) {
  if (!node) throw new Error(`Missing mount point #${name} in index.html`);
}

const app = new App({ api: timelineApi, roots });
await app.start();

// Handy in the browser console while learning: `historyTimeline.state`
globalThis.historyTimeline = app;
