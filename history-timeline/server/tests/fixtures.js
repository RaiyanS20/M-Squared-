import { Country, EventCategory, HistoricalEvent, Year } from '../src/domain/index.js';
import { InMemoryDataset } from '../src/repositories/memory/index.js';

/**
 * Test data builders.
 *
 * A builder with sensible defaults lets a test state only what it CARES about:
 * `anEvent({ year: 1969 })` reads as "an event in 1969", and the reader is not
 * distracted by a source URL irrelevant to the assertion.
 */
export const FRANCE = new Country({ id: 1, code: 'FRA', name: 'France', region: 'Europe' });
export const JAPAN = new Country({ id: 2, code: 'JPN', name: 'Japan', region: 'Asia' });
export const BRAZIL = new Country({ id: 3, code: 'BRA', name: 'Brazil', region: 'South America' });

let nextId = 100;

export function anEvent(overrides = {}) {
  const { year = 1969, ...rest } = overrides;
  return new HistoricalEvent({
    id: nextId++,
    countryId: FRANCE.id,
    title: 'A documented event',
    summary: 'Something verifiable happened.',
    category: EventCategory.Politics,
    sourceUrl: 'https://en.wikipedia.org/wiki/History',
    monthDay: null,
    ...rest,
    // `year` is rebuilt last so callers can pass a plain number.
    year: year instanceof Year ? year : new Year(year),
  });
}

/** A small world: France has two 1969 events, Japan one, Brazil none. */
export function aDataset() {
  return new InMemoryDataset(
    [FRANCE, JAPAN, BRAZIL],
    [
      anEvent({ id: 1, countryId: FRANCE.id, year: 1969, title: 'Concorde first flight', monthDay: '03-02' }),
      anEvent({ id: 2, countryId: FRANCE.id, year: 1969, title: 'De Gaulle resigns', monthDay: '04-28' }),
      anEvent({ id: 3, countryId: FRANCE.id, year: 1981, title: 'TGV service opens' }),
      anEvent({ id: 4, countryId: JAPAN.id, year: 1969, title: 'A Japanese event' }),
    ],
  );
}
