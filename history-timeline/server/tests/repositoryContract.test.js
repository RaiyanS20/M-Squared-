import { describe, it, after } from 'node:test';
import assert from 'node:assert/strict';
import {
  InMemoryCountryRepository,
  InMemoryEventRepository,
} from '../src/repositories/memory/index.js';
import {
  PostgresCountryRepository,
  PostgresEventRepository,
} from '../src/repositories/postgres/index.js';
import { PostgresDatabase } from '../src/db/Database.js';
import { aDataset } from './fixtures.js';

/**
 * A CONTRACT TEST.
 *
 * `CountryRepository` is an abstraction with two implementations. An abstraction
 * is only trustworthy if every implementation behaves the SAME WAY — otherwise
 * "it passed in tests" tells you nothing about production.
 *
 * This matters even more in JavaScript than in a typed language. A compiler
 * would at least tell you a method was missing; here, nothing checks that
 * `InMemoryCountryRepository` and `PostgresCountryRepository` agree until
 * something calls them. This suite is that something.
 *
 * So the assertions live in ONE function and each implementation is fed through
 * it. The in-memory pair always runs. The Postgres pair runs only when
 * TEST_DATABASE_URL points at a migrated, seeded database:
 *
 *   npm run db:up && npm run db:migrate && npm run db:seed
 *   TEST_DATABASE_URL=postgres://historian:historian@localhost:5432/history_timeline npm test
 *
 * Fast feedback by default; real confidence on demand. CI sets the variable, so
 * nothing is permanently skipped.
 */
function contractFor(label, build, { knownCode, populatedYear }) {
  describe(`${label} satisfies the repository contract`, () => {
    it('findAll returns countries sorted by name', async () => {
      const { countries } = await build();
      const names = (await countries.findAll()).map((c) => c.name);
      assert.deepEqual(names, [...names].sort((a, b) => a.localeCompare(b)));
      assert.ok(names.length > 0);
    });

    it('findByCode is case-insensitive', async () => {
      const { countries } = await build();
      const upper = await countries.findByCode(knownCode.toUpperCase());
      const lower = await countries.findByCode(knownCode.toLowerCase());
      assert.ok(upper, 'expected the known country to exist');
      assert.equal(lower?.id, upper.id);
    });

    it('findByCode returns null for an unknown code', async () => {
      const { countries } = await build();
      assert.equal(await countries.findByCode('ZZZ'), null);
    });

    it('findById returns null for an unknown id', async () => {
      const { countries } = await build();
      assert.equal(await countries.findById(-1), null);
    });

    it('findAllWithEventsInYear returns only countries that have events', async () => {
      const { countries, events } = await build();
      const withEvents = await countries.findAllWithEventsInYear(populatedYear);
      assert.ok(withEvents.length > 0);
      for (const country of withEvents) {
        const found = await events.findByYearAndCountry(populatedYear, country.id);
        assert.ok(found.length > 0, `${country.name} was listed but has no events`);
      }
    });

    it('findAllWithEventsInYear is empty for a year with nothing recorded', async () => {
      const { countries } = await build();
      // 1902 is deliberately absent from both the fixture and the seed dataset.
      assert.deepEqual(await countries.findAllWithEventsInYear(1902), []);
    });

    it('findByYearAndCountry returns events in date order', async () => {
      const { countries, events } = await build();
      const country = await countries.findByCode(knownCode);
      const found = await events.findByYearAndCountry(populatedYear, country.id);

      const dated = found.filter((e) => e.monthDay !== null).map((e) => e.monthDay);
      assert.deepEqual(dated, [...dated].sort());

      // Undated events must come last.
      const firstUndated = found.findIndex((e) => e.monthDay === null);
      if (firstUndated !== -1) {
        assert.ok(found.slice(firstUndated).every((e) => e.monthDay === null));
      }
    });

    it('countByYearRange returns ascending years with positive counts', async () => {
      const { events } = await build();
      const counts = await events.countByYearRange(1900, 2000);
      const years = counts.map((c) => c.year);
      assert.deepEqual(years, [...years].sort((a, b) => a - b));
      assert.ok(counts.every((c) => Number.isInteger(c.eventCount) && c.eventCount > 0));
    });

    it('countByYearRange excludes years outside the range', async () => {
      const { events } = await build();
      const counts = await events.countByYearRange(1969, 1969);
      assert.ok(counts.every((c) => c.year === 1969));
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
let sharedDb = null;

after(async () => {
  await sharedDb?.close();
});

describe('PostgreSQL repositories', { skip: !testDatabaseUrl && 'set TEST_DATABASE_URL to run' }, () => {
  contractFor(
    'Postgres repositories',
    async () => {
      sharedDb ??= new PostgresDatabase({ connectionString: testDatabaseUrl });
      return {
        countries: new PostgresCountryRepository(sharedDb),
        events: new PostgresEventRepository(sharedDb),
      };
    },
    // Matches the seed dataset: France has the Treaty of Versailles in 1919.
    { knownCode: 'FRA', populatedYear: 1919 },
  );
});
