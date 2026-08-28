import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import { Country, HistoricalEvent, Year } from '../src/domain/index.js';
import { SEED_COUNTRIES, SEED_EVENTS } from '../src/db/seedData.js';

/**
 * The CONTENT is tested, not just the code.
 *
 * For an educational product the data IS the product, so its editorial rules
 * deserve the same protection as the business rules. This suite fails the build
 * if someone adds an event with no citation, a duplicate, or an unknown country
 * code — before it ever reaches a learner.
 *
 * Whatever your product's equivalent is, test it. The rules that define
 * correctness are not always in the code.
 */
describe('seed dataset', () => {
  it('every country is valid', () => {
    for (const c of SEED_COUNTRIES) {
      assert.doesNotThrow(() => new Country({ id: 0, ...c }), `${c.name} is invalid`);
    }
  });

  it('country codes are unique', () => {
    const codes = SEED_COUNTRIES.map((c) => c.code);
    assert.equal(new Set(codes).size, codes.length);
  });

  it('every event is valid, cited, and in range', () => {
    for (const e of SEED_EVENTS) {
      assert.doesNotThrow(
        () => new HistoricalEvent({ id: 0, countryId: 0, ...e, year: new Year(e.year) }),
        `"${e.title}" (${e.countryCode} ${e.year}) is not a valid event`,
      );
    }
  });

  it('every event belongs to a known country', () => {
    const codes = new Set(SEED_COUNTRIES.map((c) => c.code));
    for (const e of SEED_EVENTS) {
      assert.ok(codes.has(e.countryCode), `unknown country code in "${e.title}"`);
    }
  });

  it('no country has a duplicate event in the same year', () => {
    const seen = new Set();
    for (const e of SEED_EVENTS) {
      const key = `${e.countryCode}|${e.year}|${e.title}`;
      assert.ok(!seen.has(key), `duplicate entry: ${key}`);
      seen.add(key);
    }
  });

  it('every event falls inside the MVP window', () => {
    for (const e of SEED_EVENTS) {
      assert.ok(e.year >= Year.MVP_START, `${e.title} is before ${Year.MVP_START}`);
      assert.ok(e.year <= Year.latestAllowed(), `${e.title} is in the future`);
    }
  });

  it('every country has at least three events, so no country is a dead end', () => {
    for (const c of SEED_COUNTRIES) {
      const count = SEED_EVENTS.filter((e) => e.countryCode === c.code).length;
      assert.ok(count >= 3, `${c.name} only has ${count} event(s)`);
    }
  });
});
