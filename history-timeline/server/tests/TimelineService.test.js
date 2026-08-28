import { describe, it, beforeEach } from 'node:test';
import assert from 'node:assert/strict';
import { NotFoundError, ValidationError, Year } from '../src/domain/index.js';
import {
  InMemoryCountryRepository,
  InMemoryEventRepository,
} from '../src/repositories/memory/index.js';
import { TimelineService } from '../src/services/index.js';
import { aDataset, FRANCE } from './fixtures.js';

/**
 * The payoff of dependency inversion: this exercises the REAL service, with the
 * real business rules, and needs no database, no server and no mocking library.
 * The whole file runs in single-digit milliseconds.
 */
describe('TimelineService', () => {
  let service;

  beforeEach(() => {
    const data = aDataset();
    service = new TimelineService(
      new InMemoryCountryRepository(data),
      new InMemoryEventRepository(data),
    );
  });

  describe('getTimelineScale — MVP feature 1', () => {
    it('defaults to 1900 through the current year', async () => {
      const scale = await service.getTimelineScale();
      assert.equal(scale.startYear, 1900);
      assert.equal(scale.endYear, Year.latestAllowed());
    });

    it('emits one tick per year, including years with no events', async () => {
      const scale = await service.getTimelineScale(1900, 1905);
      assert.equal(scale.ticks.length, 6);
      assert.ok(scale.ticks.every((t) => t.eventCount === 0));
    });

    it('counts events per year', async () => {
      const scale = await service.getTimelineScale(1969, 1969);
      assert.deepEqual(scale.ticks, [{ year: 1969, eventCount: 3 }]);
    });

    it('rejects a backwards range', async () => {
      await assert.rejects(() => service.getTimelineScale(2000, 1900), ValidationError);
    });

    it('rejects a year before the timeline starts', async () => {
      await assert.rejects(() => service.getTimelineScale(50, 1900), ValidationError);
    });
  });

  describe('listCountries — MVP feature 3', () => {
    it('returns every country alphabetically', async () => {
      const names = (await service.listCountries()).map((c) => c.name);
      assert.deepEqual(names, ['Brazil', 'France', 'Japan']);
    });
  });

  describe('listCountriesForYear — MVP features 2 and 3', () => {
    it('returns only countries with something recorded that year', async () => {
      const names = (await service.listCountriesForYear(1969)).map((c) => c.name);
      assert.deepEqual(names, ['France', 'Japan']);
    });

    it('returns an empty list for a quiet year', async () => {
      assert.deepEqual(await service.listCountriesForYear(1902), []);
    });
  });

  describe('getEventsForYearAndCountry — MVP feature 4', () => {
    it("returns that country's events for that year, in date order", async () => {
      const result = await service.getEventsForYearAndCountry(1969, 'FRA');
      assert.ok(result.country.equals(FRANCE));
      assert.deepEqual(result.events.map((e) => e.title), [
        'Concorde first flight', // 2 March
        'De Gaulle resigns', //    28 April
      ]);
    });

    it('accepts a lowercase country code', async () => {
      const result = await service.getEventsForYearAndCountry(1969, 'fra');
      assert.equal(result.events.length, 2);
    });

    it('returns an empty list — not an error — for a year with no events', async () => {
      const result = await service.getEventsForYearAndCountry(1955, 'FRA');
      assert.deepEqual(result.events, []);
    });

    it('raises NotFoundError for an unknown country', async () => {
      await assert.rejects(() => service.getEventsForYearAndCountry(1969, 'ZZZ'), NotFoundError);
    });

    it('raises ValidationError for an impossible year', async () => {
      await assert.rejects(() => service.getEventsForYearAndCountry(3000, 'FRA'), ValidationError);
    });
  });
});
