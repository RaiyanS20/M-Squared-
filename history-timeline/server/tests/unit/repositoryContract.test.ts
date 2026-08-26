import { afterAll, describe, expect, it } from 'vitest';
import { CountryRepository, EventRepository } from '../../src/repositories/index.js';
import {
  InMemoryCountryRepository,
  InMemoryEventRepository,
} from '../../src/repositories/memory/index.js';
import {
  PostgresCountryRepository,
  PostgresEventRepository,
} from '../../src/repositories/postgres/index.js';
import { PostgresDatabase } from '../../src/db/Database.js';
import { aDataset } from '../fixtures.js';

/**
 * A CONTRACT TEST.
 *
 * `CountryRepository` is an abstraction with two implementations. An abstraction
 * is only trustworthy if every implementation behaves the same way — otherwise
 * "it passed in tests" means nothing about production.
 *
 * So the assertions live in one function, and each implementation is fed through
 * it. The in-memory pair always runs. The Postgres pair runs only when
 * TEST_DATABASE_URL points at a migrated, seeded database:
 *
 *     npm run db:up && npm run db:migrate && npm run db:seed
 *     TEST_DATABASE_URL=postgres://historian:historian@localhost:5432/history_timeline npm test
 *
 * This is how you get fast feedback by default AND real confidence on demand.
 */
function contractFor(
  label: string,
  build: () => Promise<{ countries: CountryRepository; events: EventRepository }>,
  expectations: {
    /** A country code known to exist in this fixture. */
    knownCode: string;
    /** A year in which `knownCode` has at least one event. */
    populatedYear: number;
  },
) {
  describe(`${label} satisfies the repository contract`, () => {
    it('findAll returns countries sorted by name', async () => {
      const { countries } = await build();
      const names = (await countries.findAll()).map((c) => c.name);
      expect(names).toEqual([...names].sort((a, b) => a.localeCompare(b)));
      expect(names.length).toBeGreaterThan(0);
    });

    it('findByCode is case-insensitive', async () => {
      const { countries } = await build();
      const upper = await countries.findByCode(expectations.knownCode.toUpperCase());
      const lower = await countries.findByCode(expectations.knownCode.toLowerCase());
      expect(upper).not.toBeNull();
      expect(lower?.id).toBe(upper?.id);
    });

    it('findByCode returns null for an unknown code', async () => {
      const { countries } = await build();
      expect(await countries.findByCode('ZZZ')).toBeNull();
    });

    it('findById returns null for an unknown id', async () => {
      const { countries } = await build();
      expect(await countries.findById(-1)).toBeNull();
    });

    it('findAllWithEventsInYear returns only countries that have events', async () => {
      const { countries, events } = await build();
      const withEvents = await countries.findAllWithEventsInYear(expectations.populatedYear);
      expect(withEvents.length).toBeGreaterThan(0);
      for (const country of withEvents) {
        const found = await events.findByYearAndCountry(expectations.populatedYear, country.id);
        expect(found.length).toBeGreaterThan(0);
      }
    });

    it('findAllWithEventsInYear is empty for a year with nothing recorded', async () => {
      const { countries } = await build();
      // 1902 is deliberately absent from both the test fixture and the seed dataset.
      expect(await countries.findAllWithEventsInYear(1902)).toEqual([]);
    });

    it('findByYearAndCountry returns events in date order', async () => {
      const { countries, events } = await build();
      const country = await countries.findByCode(expectations.knownCode);
      const found = await events.findByYearAndCountry(expectations.populatedYear, country!.id);
      const dated = found.filter((e) => e.monthDay !== null).map((e) => e.monthDay!);
      expect(dated).toEqual([...dated].sort());
      // Undated events must come last.
      const firstUndated = found.findIndex((e) => e.monthDay === null);
      if (firstUndated !== -1) {
        expect(found.slice(firstUndated).every((e) => e.monthDay === null)).toBe(true);
      }
    });

    it('countByYearRange returns ascending years with positive counts', async () => {
      const { events } = await build();
      const counts = await events.countByYearRange(1900, 2000);
      expect(counts.map((c) => c.year)).toEqual([...counts.map((c) => c.year)].sort((a, b) => a - b));
      expect(counts.every((c) => Number.isInteger(c.eventCount) && c.eventCount > 0)).toBe(true);
    });

    it('countByYearRange excludes years outside the range', async () => {
      const { events } = await build();
      const counts = await events.countByYearRange(1969, 1969);
      expect(counts.every((c) => c.year === 1969)).toBe(true);
    });
  });
}

// ---- Implementation 1: in-memory (always runs) ------------------------------
contractFor(
  'InMemory repositories',
  async () => {
    const data = aDataset();
    return {
      countries: new InMemoryCountryRepository(data),
      events: new InMemoryEventRepository(data),
    };
  },
  { knownCode: 'FRA', populatedYear: 1969 },
);

// ---- Implementation 2: PostgreSQL (opt-in) ----------------------------------
const testDatabaseUrl = process.env.TEST_DATABASE_URL;
let sharedDb: PostgresDatabase | null = null;

afterAll(async () => {
  await sharedDb?.close();
});

describe.skipIf(!testDatabaseUrl)('PostgreSQL repositories', () => {
  contractFor(
    'Postgres repositories',
    async () => {
      sharedDb ??= new PostgresDatabase({ connectionString: testDatabaseUrl! });
      return {
        countries: new PostgresCountryRepository(sharedDb),
        events: new PostgresEventRepository(sharedDb),
      };
    },
    // Matches the seed dataset: France has two events in 1969... adjust if the
    // seed changes. 1919 has the Treaty of Versailles for FRA.
    { knownCode: 'FRA', populatedYear: 1919 },
  );
});
