import { describe, expect, it } from 'vitest';
import { Country, EventCategory, HistoricalEvent, ValidationError, Year } from '../../src/domain/index.js';
import { anEvent } from '../fixtures.js';

describe('Country', () => {
  it('normalises the code to uppercase and trims the name', () => {
    const c = new Country({ id: 1, code: 'fra', name: '  France ', region: 'Europe' });
    expect(c.code).toBe('FRA');
    expect(c.name).toBe('France');
  });

  it.each(['FR', 'FRAN', '12A', ''])('rejects the invalid code %o', (code) => {
    expect(() => new Country({ id: 1, code, name: 'X', region: 'Europe' })).toThrow(ValidationError);
  });

  it('rejects a blank name', () => {
    expect(() => new Country({ id: 1, code: 'FRA', name: '   ', region: 'Europe' })).toThrow(
      ValidationError,
    );
  });

  it('defaults an unknown region rather than failing', () => {
    expect(new Country({ id: 1, code: 'FRA', name: 'France', region: '' }).region).toBe('Unknown');
  });

  // Entities compare by identity: this is the defining difference from Year.
  it('treats two objects with the same id as the same country', () => {
    const a = new Country({ id: 7, code: 'FRA', name: 'France', region: 'Europe' });
    const b = new Country({ id: 7, code: 'FRA', name: 'French Republic', region: 'Europe' });
    expect(a.equals(b)).toBe(true);
  });
});

describe('HistoricalEvent', () => {
  it('builds a valid event', () => {
    const e = anEvent({ year: 1969, title: 'Apollo 11' });
    expect(e.year.value).toBe(1969);
    expect(e.title).toBe('Apollo 11');
  });

  // The editorial rule of the whole project, enforced in code.
  it.each(['', 'not-a-url', 'ftp://example.com/x', 'wikipedia.org/wiki/X'])(
    'refuses to exist without a citable source (%o)',
    (sourceUrl) => {
      expect(() => anEvent({ sourceUrl })).toThrow(ValidationError);
    },
  );

  it('rejects a title longer than the limit', () => {
    expect(() => anEvent({ title: 'x'.repeat(HistoricalEvent.MAX_TITLE_LENGTH + 1) })).toThrow(
      ValidationError,
    );
  });

  it('rejects a summary longer than the limit', () => {
    expect(() => anEvent({ summary: 'x'.repeat(HistoricalEvent.MAX_SUMMARY_LENGTH + 1) })).toThrow(
      ValidationError,
    );
  });

  it('rejects an unknown category', () => {
    expect(() => anEvent({ category: 'gossip' as EventCategory })).toThrow(ValidationError);
  });

  it.each(['3-2', '0302', '03/02'])('rejects the malformed monthDay %o', (monthDay) => {
    expect(() => anEvent({ monthDay })).toThrow(ValidationError);
  });

  it('allows a year-only event', () => {
    expect(anEvent({ monthDay: null }).monthDay).toBeNull();
  });

  describe('compareByDate', () => {
    it('orders by year first', () => {
      expect(anEvent({ year: 1900 }).compareByDate(anEvent({ year: 1901 }))).toBeLessThan(0);
    });

    it('orders by month and day within a year', () => {
      const march = anEvent({ year: 1969, monthDay: '03-02' });
      const april = anEvent({ year: 1969, monthDay: '04-28' });
      expect(march.compareByDate(april)).toBeLessThan(0);
    });

    it('puts dated events before undated ones', () => {
      const dated = anEvent({ year: 1969, monthDay: '03-02' });
      const undated = anEvent({ year: 1969, monthDay: null });
      expect(dated.compareByDate(undated)).toBeLessThan(0);
      expect(undated.compareByDate(dated)).toBeGreaterThan(0);
    });

    it('falls back to title for two undated events', () => {
      const a = anEvent({ year: 1969, monthDay: null, title: 'Alpha' });
      const b = anEvent({ year: 1969, monthDay: null, title: 'Beta' });
      expect(a.compareByDate(b)).toBeLessThan(0);
    });
  });

  it('serialises the year as a number for the API', () => {
    const json = anEvent({ year: 1969 }).toJSON();
    expect(json.year).toBe(1969);
    expect(new Year(json.year).value).toBe(1969);
  });
});
