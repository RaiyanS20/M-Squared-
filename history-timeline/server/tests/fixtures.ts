import { Country, EventCategory, HistoricalEvent, Year } from '../src/domain/index.js';
import { InMemoryDataset } from '../src/repositories/memory/index.js';

/**
 * Test data builders.
 *
 * A builder with sensible defaults means a test only states what it CARES about:
 * `anEvent({ year: 1969 })` reads as "an event in 1969", and the reader is not
 * distracted by a source URL that is irrelevant to the assertion.
 */
export const FRANCE = new Country({ id: 1, code: 'FRA', name: 'France', region: 'Europe' });
export const JAPAN = new Country({ id: 2, code: 'JPN', name: 'Japan', region: 'Asia' });
export const BRAZIL = new Country({ id: 3, code: 'BRA', name: 'Brazil', region: 'South America' });

let nextId = 100;

export function anEvent(overrides: Partial<{
  id: number;
  countryId: number;
  year: number;
  title: string;
  summary: string;
  category: EventCategory;
  sourceUrl: string;
  monthDay: string | null;
}> = {}): HistoricalEvent {
  return new HistoricalEvent({
    id: overrides.id ?? nextId++,
    countryId: overrides.countryId ?? FRANCE.id,
    year: new Year(overrides.year ?? 1969),
    title: overrides.title ?? 'A documented event',
    summary: overrides.summary ?? 'Something verifiable happened.',
    category: overrides.category ?? EventCategory.Politics,
    sourceUrl: overrides.sourceUrl ?? 'https://en.wikipedia.org/wiki/History',
    monthDay: overrides.monthDay === undefined ? null : overrides.monthDay,
  });
}

/** A small world: France has two 1969 events, Japan one, Brazil none. */
export function aDataset(): InMemoryDataset {
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
