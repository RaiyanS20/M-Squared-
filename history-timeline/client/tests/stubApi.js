/** Shared browser-side test data, mirroring what the API returns. */
export const FRANCE = { id: 1, code: 'FRA', name: 'France', region: 'Europe' };
export const JAPAN = { id: 2, code: 'JPN', name: 'Japan', region: 'Asia' };

export function anEvent(overrides = {}) {
  return {
    id: 1,
    countryId: FRANCE.id,
    year: 1969,
    title: 'Concorde first flight',
    summary: 'The prototype flew from Toulouse.',
    category: 'science',
    sourceUrl: 'https://en.wikipedia.org/wiki/Concorde',
    monthDay: '03-02',
    ...overrides,
  };
}

export const TIMELINE = {
  startYear: 1967,
  endYear: 1970,
  ticks: [
    { year: 1967, eventCount: 0 },
    { year: 1968, eventCount: 1 },
    { year: 1969, eventCount: 2 },
    { year: 1970, eventCount: 0 },
  ],
};

export const COUNTRIES_1969 = { year: 1969, countries: [FRANCE, JAPAN] };

export const EVENTS_FRA_1969 = {
  year: 1969,
  country: FRANCE,
  events: [
    anEvent(),
    anEvent({ id: 2, title: 'De Gaulle resigns', category: 'politics', monthDay: '04-28' }),
  ],
};

/**
 * Replaces `globalThis.fetch` with a tiny router, and returns a record of the
 * URLs that were requested.
 *
 * We stub at the FETCH boundary rather than replacing `TimelineApi`, so
 * `ApiClient`'s real error handling, JSON parsing and abort logic all run during
 * the tests. Mocking our own class would mean the tests never execute the code
 * we actually ship.
 */
export function stubFetch(routes, { status = {} } = {}) {
  const calls = [];
  const original = globalThis.fetch;

  globalThis.fetch = async (input, init = {}) => {
    const url = String(input);
    calls.push(url);

    // Honour an abort signal, so tests can exercise cancellation for real.
    if (init.signal?.aborted) {
      const error = new Error('The operation was aborted.');
      error.name = 'AbortError';
      throw error;
    }

    const match = Object.keys(routes).find((route) => url.includes(route));
    if (!match) {
      return new Response(
        JSON.stringify({ error: { code: 'ROUTE_NOT_FOUND', message: `No stub for ${url}` } }),
        { status: 404, headers: { 'Content-Type': 'application/json' } },
      );
    }

    const body = routes[match];
    return new Response(typeof body === 'string' ? body : JSON.stringify(body), {
      status: status[match] ?? 200,
      headers: { 'Content-Type': 'application/json' },
    });
  };

  return { calls, restore: () => (globalThis.fetch = original) };
}
