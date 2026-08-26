import { vi } from 'vitest';
import type {
  CountriesForYear,
  Country,
  EventsForYearAndCountry,
  HistoricalEvent,
  TimelineScale,
} from '../src/api/index.js';

export const FRANCE: Country = { id: 1, code: 'FRA', name: 'France', region: 'Europe' };
export const JAPAN: Country = { id: 2, code: 'JPN', name: 'Japan', region: 'Asia' };

export function anEvent(overrides: Partial<HistoricalEvent> = {}): HistoricalEvent {
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

export const TIMELINE: TimelineScale = {
  startYear: 1967,
  endYear: 1970,
  ticks: [
    { year: 1967, eventCount: 0 },
    { year: 1968, eventCount: 1 },
    { year: 1969, eventCount: 2 },
    { year: 1970, eventCount: 0 },
  ],
};

export const COUNTRIES_1969: CountriesForYear = { year: 1969, countries: [FRANCE, JAPAN] };

export const EVENTS_FRA_1969: EventsForYearAndCountry = {
  year: 1969,
  country: FRANCE,
  events: [anEvent(), anEvent({ id: 2, title: 'De Gaulle resigns', category: 'politics', monthDay: '04-28' })],
};

/**
 * Stubs `globalThis.fetch` with a tiny router.
 *
 * We stub at the FETCH boundary rather than mocking our own `TimelineApi`. That
 * means `ApiClient`'s real error handling, JSON parsing and abort logic all run
 * during these tests — so the tests cover the code we actually ship, and would
 * catch a bug in the client's own HTTP layer.
 */
export function stubFetch(
  routes: Record<string, unknown>,
  options: { status?: Record<string, number> } = {},
) {
  return vi.spyOn(globalThis, 'fetch').mockImplementation(async (input) => {
    const url = typeof input === 'string' ? input : (input as Request).url ?? String(input);
    const match = Object.keys(routes).find((route) => url.includes(route));

    if (!match) {
      return new Response(
        JSON.stringify({ error: { code: 'ROUTE_NOT_FOUND', message: `No stub for ${url}` } }),
        { status: 404, headers: { 'Content-Type': 'application/json' } },
      );
    }

    const status = options.status?.[match] ?? 200;
    return new Response(JSON.stringify(routes[match]), {
      status,
      headers: { 'Content-Type': 'application/json' },
    });
  });
}
