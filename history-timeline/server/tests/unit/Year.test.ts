import { describe, expect, it } from 'vitest';
import { ValidationError, Year } from '../../src/domain/index.js';

/**
 * Testing lesson: test the RULES, not the implementation.
 * Every assertion here corresponds to a sentence you could say out loud about
 * the product ("the timeline cannot show future years"). If the internals of
 * Year are rewritten, these tests should still pass unchanged.
 */
describe('Year', () => {
  describe('construction', () => {
    it('accepts a year inside the supported span', () => {
      expect(new Year(1969).value).toBe(1969);
    });

    it('accepts the earliest supported year, 100 AD', () => {
      expect(new Year(Year.EARLIEST).value).toBe(100);
    });

    it('rejects years before 100 AD', () => {
      expect(() => new Year(99)).toThrow(ValidationError);
    });

    it('rejects future years', () => {
      expect(() => new Year(Year.latestAllowed() + 1)).toThrow(ValidationError);
    });

    it('rejects non-integers', () => {
      expect(() => new Year(1969.5)).toThrow(ValidationError);
      expect(() => new Year(Number.NaN)).toThrow(ValidationError);
    });

    // Boundary values are where bugs live. Always test the edge, the edge minus
    // one, and the edge plus one.
    it('accepts the current year', () => {
      expect(() => new Year(Year.latestAllowed())).not.toThrow();
    });
  });

  describe('fromString', () => {
    it('parses a numeric string', () => {
      expect(Year.fromString(' 1900 ').value).toBe(1900);
    });

    it.each(['', 'nineteen-sixty-nine', '19a9', '1969.0'])(
      'rejects %o',
      (raw) => {
        expect(() => Year.fromString(raw)).toThrow(ValidationError);
      },
    );
  });

  describe('derived values', () => {
    it('reports the decade', () => {
      expect(new Year(1917).decade).toBe(1910);
      expect(new Year(1900).decade).toBe(1900);
    });

    it('reports the century, 1-based', () => {
      expect(new Year(1917).century).toBe(20);
      expect(new Year(1900).century).toBe(19);
      expect(new Year(2001).century).toBe(21);
    });
  });

  describe('comparison', () => {
    it('compares by value, not identity', () => {
      expect(new Year(1969).equals(new Year(1969))).toBe(true);
      expect(new Year(1969) === new Year(1969)).toBe(false);
    });

    it('orders years', () => {
      expect(new Year(1900).isBefore(new Year(1901))).toBe(true);
      expect(new Year(1901).isAfter(new Year(1900))).toBe(true);
    });
  });

  describe('range', () => {
    it('is inclusive at both ends', () => {
      const years = Year.range(new Year(1900), new Year(1903));
      expect(years.map((y) => y.value)).toEqual([1900, 1901, 1902, 1903]);
    });

    it('returns a single year when start equals end', () => {
      expect(Year.range(new Year(1969), new Year(1969))).toHaveLength(1);
    });

    it('refuses a backwards range', () => {
      expect(() => Year.range(new Year(2000), new Year(1900))).toThrow(ValidationError);
    });

    it('mvpRange starts at 1900 and ends this year', () => {
      const range = Year.mvpRange();
      expect(range[0]!.value).toBe(1900);
      expect(range.at(-1)!.value).toBe(Year.latestAllowed());
    });
  });

  it('serialises to a plain number', () => {
    expect(JSON.parse(JSON.stringify({ y: new Year(1969) }))).toEqual({ y: 1969 });
  });

  it('is immutable', () => {
    const year = new Year(1969);
    expect(() => {
      (year as unknown as { value: number }).value = 2000;
    }).toThrow();
  });
});
