import { describe, it, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { createApp } from '../src/api/index.js';
import {
  InMemoryCountryRepository,
  InMemoryEventRepository,
} from '../src/repositories/memory/index.js';
import { TimelineService } from '../src/services/index.js';
import { Year } from '../src/domain/index.js';
import { aDataset } from './fixtures.js';

/**
 * INTEGRATION TESTS: the full HTTP stack — routing, validation, controllers,
 * services, repositories, error handling — with in-memory data underneath.
 *
 * This layer catches what unit tests cannot see: a route registered at the wrong
 * path, a 500 where a 404 belongs, a field renamed in the JSON that the browser
 * app still expects.
 *
 * No supertest, no test-HTTP library. We start the real app on port 0 — "any
 * free port" — and talk to it with the `fetch` that is built into Node. Fewer
 * dependencies, and it makes clear that an HTTP test is just an HTTP request.
 */
let server;
let baseUrl;

before(async () => {
  const data = aDataset();
  const service = new TimelineService(
    new InMemoryCountryRepository(data),
    new InMemoryEventRepository(data),
  );
  // serveClient: false keeps the test focused on the API. With it on, an unknown
  // path would fall through to the static handler instead of the JSON 404.
  const app = createApp({ service, db: null, isProduction: false, serveClient: false });

  await new Promise((resolve) => {
    server = app.listen(0, resolve);
  });
  baseUrl = `http://127.0.0.1:${server.address().port}`;
});

after(async () => {
  await new Promise((resolve) => server.close(resolve));
});

/** Small helper so each test reads as one line of intent. */
async function get(path) {
  const response = await fetch(`${baseUrl}${path}`);
  return { status: response.status, body: await response.json(), headers: response.headers };
}

describe('GET /health and /ready', () => {
  it('reports ok', async () => {
    const { status, body } = await get('/health');
    assert.equal(status, 200);
    assert.equal(body.status, 'ok');
  });

  it('reports ready when no database is configured', async () => {
    const { status, body } = await get('/ready');
    assert.equal(status, 200);
    assert.equal(body.database, 'not-configured');
  });
});

describe('GET /api/timeline', () => {
  it('returns the default MVP scale', async () => {
    const { status, body } = await get('/api/timeline');
    assert.equal(status, 200);
    assert.equal(body.startYear, 1900);
    assert.equal(body.endYear, Year.latestAllowed());
    assert.deepEqual(body.ticks[0], { year: 1900, eventCount: 0 });
  });

  it('honours an explicit range', async () => {
    const { body } = await get('/api/timeline?startYear=1969&endYear=1970');
    assert.deepEqual(body.ticks, [
      { year: 1969, eventCount: 3 },
      { year: 1970, eventCount: 0 },
    ]);
  });

  it('rejects a non-numeric year with 400 and a machine-readable code', async () => {
    const { status, body } = await get('/api/timeline?startYear=abc');
    assert.equal(status, 400);
    assert.equal(body.error.code, 'VALIDATION_ERROR');
    assert.match(body.error.message, /not a valid year/i);
  });

  it('rejects a backwards range with 400', async () => {
    assert.equal((await get('/api/timeline?startYear=2000&endYear=1900')).status, 400);
  });

  it('rejects a repeated query parameter rather than misbehaving', async () => {
    assert.equal((await get('/api/timeline?startYear=1900&startYear=1950')).status, 400);
  });
});

describe('GET /api/countries', () => {
  it('lists every country', async () => {
    const { body } = await get('/api/countries');
    assert.deepEqual(body.countries.map((c) => c.code), ['BRA', 'FRA', 'JPN']);
  });

  // Pins the public contract: this fails if a field is ever added or removed,
  // which is exactly what you want from an API other programs depend on.
  it('exposes only the intended fields', async () => {
    const { body } = await get('/api/countries');
    assert.deepEqual(Object.keys(body.countries[0]).sort(), ['code', 'id', 'name', 'region']);
  });
});

describe('GET /api/years/:year/countries', () => {
  it('lists only countries with events in that year', async () => {
    const { body } = await get('/api/years/1969/countries');
    assert.equal(body.year, 1969);
    assert.deepEqual(body.countries.map((c) => c.code), ['FRA', 'JPN']);
  });

  it('returns an empty list for a quiet year', async () => {
    const { status, body } = await get('/api/years/1902/countries');
    assert.equal(status, 200);
    assert.deepEqual(body.countries, []);
  });

  it('rejects a year before 100 AD', async () => {
    assert.equal((await get('/api/years/50/countries')).status, 400);
  });

  it('rejects a future year', async () => {
    assert.equal((await get(`/api/years/${Year.latestAllowed() + 1}/countries`)).status, 400);
  });
});

describe('GET /api/years/:year/countries/:code/events', () => {
  it('returns the events for that year and country, in date order', async () => {
    const { body } = await get('/api/years/1969/countries/FRA/events');
    assert.equal(body.year, 1969);
    assert.equal(body.country.name, 'France');
    assert.deepEqual(body.events.map((e) => e.title), [
      'Concorde first flight',
      'De Gaulle resigns',
    ]);
  });

  it('serialises an event with the fields the browser app needs', async () => {
    const { body } = await get('/api/years/1969/countries/FRA/events');
    assert.deepEqual(Object.keys(body.events[0]).sort(), [
      'category', 'countryId', 'id', 'monthDay', 'sourceUrl', 'summary', 'title', 'year',
    ]);
  });

  it('accepts a lowercase code', async () => {
    assert.equal((await get('/api/years/1969/countries/fra/events')).status, 200);
  });

  it('returns 404 for an unknown country', async () => {
    const { status, body } = await get('/api/years/1969/countries/ZZZ/events');
    assert.equal(status, 404);
    assert.equal(body.error.code, 'NOT_FOUND');
  });

  it('returns 400 for a malformed country code', async () => {
    const { status, body } = await get('/api/years/1969/countries/FRANCE/events');
    assert.equal(status, 400);
    assert.equal(body.error.code, 'VALIDATION_ERROR');
  });

  // The distinction that confuses people: "nothing recorded" is a valid answer
  // to a valid question, so it is 200 with an empty list. 404 would mean France
  // does not exist.
  it('returns 200 with an empty list when nothing was recorded', async () => {
    const { status, body } = await get('/api/years/1955/countries/FRA/events');
    assert.equal(status, 200);
    assert.deepEqual(body.events, []);
  });
});

describe('unknown routes and headers', () => {
  it('return JSON, not an HTML error page', async () => {
    const { status, body } = await get('/api/nope');
    assert.equal(status, 404);
    assert.equal(body.error.code, 'ROUTE_NOT_FOUND');
  });

  it('does not advertise the framework', async () => {
    const { headers } = await get('/health');
    assert.equal(headers.get('x-powered-by'), null);
  });
});
