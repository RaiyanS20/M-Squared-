import { beforeEach, describe, expect, it } from 'vitest';
import { NotFoundError, ValidationError, Year } from '../../src/domain/index.js';
import {
  InMemoryCountryRepository,
  InMemoryEventRepository,
} from '../../src/repositories/memory/index.js';
import { TimelineService } from '../../src/services/index.js';
import { aDataset, FRANCE } from '../fixtures.js';

/**
 * The payoff of dependency inversion: this exercises the real service, with the
 * real business rules, and needs no database, no server and no mocking library.
 * The whole file runs in single-digit milliseconds.
 */
describe('TimelineService', () => {
  let service: TimelineService;

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
      expect(scale.startYear).toBe(1900);
      expect(scale.endYear).toBe(Year.latestAllowed());
    });

    it('emits one tick per year, including years with no events', async () => {
      const scale = await service.getTimelineScale(1900, 1905);
      expect(scale.ticks).toHaveLength(6);
      expect(scale.ticks.every((t) => t.eventCount === 0)).toBe(true);
    });

    it('counts events per year', async () => {
      const scale = await service.getTimelineScale(1969, 1969);
      expect(scale.ticks).toEqual([{ year: 1969, eventCount: 3 }]);
    });

    it('rejects a backwards range', async () => {
      await expect(service.getTimelineScale(2000, 1900)).rejects.toThrow(ValidationError);
    });

    it('rejects a year before the timeline starts', async () => {
      await expect(service.getTimelineScale(50, 1900)).rejects.toThrow(ValidationError);
    });
  });

  describe('listCountries — MVP feature 3', () => {
    it('returns every country alphabetically', async () => {
      const names = (await service.listCountries()).map((c) => c.name);
      expect(names).toEqual(['Brazil', 'France', 'Japan']);
    });
  });

  describe('listCountriesForYear — MVP features 2 and 3', () => {
    it('returns only countries with something recorded that year', async () => {
      const names = (await service.listCountriesForYear(1969)).map((c) => c.name);
      expect(names).toEqual(['France', 'Japan']);
      expect(names).not.toContain('Brazil');
    });

    it('returns an empty list for a year with no events', async () => {
      expect(await service.listCountriesForYear(1904)).toEqual([]);
    });
  });

  describe('getEventsForYearAndCountry — MVP feature 4', () => {
    it('returns that country\'s events for that year, in date order', async () => {
      const result = await service.getEventsForYearAndCountry(1969, 'FRA');
      expect(result.country.equals(FRANCE)).toBe(true);
      expect(result.events.map((e) => e.title)).toEqual([
        'Concorde first flight', // 02 March
        'De Gaulle resigns', //    28 April
      ]);
    });

    it('accepts a lowercase country code', async () => {
      const result = await service.getEventsForYearAndCountry(1969, 'fra');
      expect(result.events).toHaveLength(2);
    });

    it('returns an empty list — not an error — for a year with no events', async () => {
      const result = await service.getEventsForYearAndCountry(1955, 'FRA');
      expect(result.events).toEqual([]);
    });

    it('raises NotFoundError for an unknown country', async () => {
      await expect(service.getEventsForYearAndCountry(1969, 'ZZZ')).rejects.toThrow(NotFoundError);
    });

    it('raises ValidationError for an impossible year', async () => {
      await expect(service.getEventsForYearAndCountry(3000, 'FRA')).rejects.toThrow(ValidationError);
    });
  });
});
