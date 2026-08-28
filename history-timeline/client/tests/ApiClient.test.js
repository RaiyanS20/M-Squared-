import { describe, it, afterEach } from 'node:test';
import assert from 'node:assert/strict';
import { ApiClient, ApiError } from '../js/api/ApiClient.js';
import { TimelineApi } from '../js/api/TimelineApi.js';
import { COUNTRIES_1969, TIMELINE, stubFetch } from './stubApi.js';

// `fetch` and `Response` are built into Node 18+, so these tests need no DOM.
let active = null;
afterEach(() => {
  active?.restore();
  active = null;
});

describe('ApiClient', () => {
  it('cannot be constructed directly — it is abstract', () => {
    assert.throws(() => new ApiClient(''), TypeError);
  });
});

describe('TimelineApi', () => {
  it('requests the timeline and returns the parsed body', async () => {
    active = stubFetch({ '/api/timeline': TIMELINE });
    assert.deepEqual(await new TimelineApi().getTimeline(), TIMELINE);
    assert.deepEqual(active.calls, ['/api/timeline']);
  });

  it('appends only the query parameters that were provided', async () => {
    active = stubFetch({ '/api/timeline': TIMELINE });
    await new TimelineApi().getTimeline(1900, undefined);
    assert.equal(active.calls[0], '/api/timeline?startYear=1900');
  });

  it('unwraps the countries array', async () => {
    active = stubFetch({ '/api/countries': { countries: COUNTRIES_1969.countries } });
    assert.equal((await new TimelineApi().getAllCountries()).length, 2);
  });

  it('honours a configured base URL', async () => {
    active = stubFetch({ '/api/timeline': TIMELINE });
    await new TimelineApi('https://api.example.com').getTimeline();
    assert.equal(active.calls[0], 'https://api.example.com/api/timeline');
  });

  it('builds the events URL from the year and country code', async () => {
    active = stubFetch({ '/api/years': { year: 1969, country: {}, events: [] } });
    await new TimelineApi().getEvents(1969, 'FRA');
    assert.equal(active.calls[0], '/api/years/1969/countries/FRA/events');
  });
});

describe('TimelineApi error handling', () => {
  it('turns a 404 body into an ApiError carrying the server code', async () => {
    active = stubFetch(
      { '/api/years': { error: { code: 'NOT_FOUND', message: "Country 'ZZZ' was not found." } } },
      { status: { '/api/years': 404 } },
    );

    const error = await new TimelineApi().getEvents(1969, 'ZZZ').catch((e) => e);
    assert.ok(error instanceof ApiError);
    assert.equal(error.code, 'NOT_FOUND');
    assert.equal(error.status, 404);
    assert.equal(error.isClientError, true);
  });

  // A crashed process or a proxy returning HTML is a real production event.
  it('survives an error response that is not JSON', async () => {
    active = stubFetch({ '/api/timeline': '<html>502 Bad Gateway</html>' }, { status: { '/api/timeline': 502 } });

    const error = await new TimelineApi().getTimeline().catch((e) => e);
    assert.ok(error instanceof ApiError);
    assert.equal(error.code, 'UNKNOWN_ERROR');
    assert.equal(error.isClientError, false); // 5xx — worth retrying
  });

  it('reports a network failure as NETWORK_ERROR', async () => {
    const original = globalThis.fetch;
    globalThis.fetch = async () => {
      throw new TypeError('Failed to fetch');
    };
    try {
      const error = await new TimelineApi().getTimeline().catch((e) => e);
      assert.equal(error.code, 'NETWORK_ERROR');
      assert.match(error.message, /could not reach/i);
    } finally {
      globalThis.fetch = original;
    }
  });

  // An abort is a normal part of the lifecycle, not a failure to report.
  it('re-throws an abort untouched so callers can ignore it', async () => {
    active = stubFetch({ '/api/timeline': TIMELINE });
    const controller = new AbortController();
    controller.abort();

    const error = await new TimelineApi().getTimeline(undefined, undefined, controller.signal).catch((e) => e);
    assert.equal(error.name, 'AbortError');
    assert.ok(!(error instanceof ApiError));
  });
});
