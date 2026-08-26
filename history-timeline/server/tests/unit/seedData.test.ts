import { describe, expect, it } from 'vitest';
import { Country, HistoricalEvent, Year } from '../../src/domain/index.js';
import { SEED_COUNTRIES, SEED_EVENTS } from '../../src/db/seedData.js';

/**
 * The content itself is tested, not just the code.
 *
 * For an educational product the data IS the product, so its editorial rules
 * deserve the same protection as the business rules. This suite fails the build
 * if someone adds an event with no citation, a duplicate entry, or a country
 * code that does not exist — before it ever reaches a learner.
 */
describe('seed dataset', () => {
  it('every country is valid', () => {
    for (const c of SEED_COUNTRIES) {
      expect(() => new Country({ id: 0, ...c })).not.toThrow();
    }
  });

  it('country codes are unique', () => {
    const codes = SEED_COUNTRIES.map((c) => c.code);
    expect(new Set(codes).size).toBe(codes.length);
  });

  it('every event is valid, cited, and in range', () => {
    for (const e of SEED_EVENTS) {
      expect(
        () =>
          new HistoricalEvent({
            id: 0,
            countryId: 0,
            year: new Year(e.year),
            title: e.title,
            summary: e.summary,
            category: e.category,
            sourceUrl: e.sourceUrl,
            monthDay: e.monthDay,
          }),
        `"${e.title}" (${e.countryCode} ${e.year}) is not a valid event`,
      ).not.toThrow();
    }
  });

  it('every event belongs to a known country', () => {
    const codes = new Set(SEED_COUNTRIES.map((c) => c.code));
    for (const e of SEED_EVENTS) {
      expect(codes, `unknown country code in "${e.title}"`).toContain(e.countryCode);
    }
  });

  it('no country has a duplicate event in the same year', () => {
    const seen = new Set<string>();
    for (const e of SEED_EVENTS) {
      const key = `${e.countryCode}|${e.year}|${e.title}`;
      expect(seen, `duplicate entry: ${key}`).not.toContain(key);
      seen.add(key);
    }
  });

  it('every event falls inside the MVP window', () => {
    for (const e of SEED_EVENTS) {
      expect(e.year).toBeGreaterThanOrEqual(Year.MVP_START);
      expect(e.year).toBeLessThanOrEqual(Year.latestAllowed());
    }
  });

  it('every country has at least three events, so no country is a dead end', () => {
    for (const c of SEED_COUNTRIES) {
      const count = SEED_EVENTS.filter((e) => e.countryCode === c.code).length;
      expect(count, `${c.name} only has ${count} event(s)`).toBeGreaterThanOrEqual(3);
    }
  });
});
