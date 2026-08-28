import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import {
  Country,
  DomainError,
  HistoricalEvent,
  NotFoundError,
  ValidationError,
  Year,
} from '../src/domain/index.js';
import { anEvent } from './fixtures.js';

/**
 * These use Node's BUILT-IN test runner — no Jest, no Vitest, no config file.
 * `node --test` finds every *.test.js and runs it. One fewer tool to learn while
 * you are learning everything else, and it makes the point that testing is not
 * magic: a test is a function that throws when something is wrong.
 *
 * Testing lesson: test the RULES, not the implementation. Every assertion here
 * corresponds to a sentence you could say out loud about the product.
 */

describe('Year', () => {
  describe('construction', () => {
    it('accepts a year inside the supported span', () => {
      assert.equal(new Year(1969).value, 1969);
    });

    // BOUNDARY TESTING: bugs live at edges, so test the edge, one below, and
    // one above. This is where off-by-one errors are caught.
    it('accepts the earliest supported year, 100 AD', () => {
      assert.equal(new Year(Year.EARLIEST).value, 100);
    });

    it('rejects years before 100 AD', () => {
      assert.throws(() => new Year(99), ValidationError);
    });

    it('accepts the current year', () => {
      assert.doesNotThrow(() => new Year(Year.latestAllowed()));
    });

    it('rejects future years', () => {
      assert.throws(() => new Year(Year.latestAllowed() + 1), ValidationError);
    });

    it('rejects values that are not whole numbers', () => {
      for (const bad of [1969.5, Number.NaN, '1969', null, undefined, {}]) {
        assert.throws(() => new Year(bad), ValidationError, `should reject ${JSON.stringify(bad)}`);
      }
    });
  });

  describe('fromString', () => {
    it('parses a numeric string, ignoring surrounding space', () => {
      assert.equal(Year.fromString(' 1900 ').value, 1900);
    });

    it('rejects anything that is not digits', () => {
      for (const bad of ['', 'nineteen-sixty-nine', '19a9', '1969.0', null]) {
        assert.throws(() => Year.fromString(bad), ValidationError, `should reject ${bad}`);
      }
    });
  });

  describe('derived values', () => {
    it('reports the decade', () => {
      assert.equal(new Year(1917).decade, 1910);
      assert.equal(new Year(1900).decade, 1900);
    });

    it('reports the century, 1-based', () => {
      assert.equal(new Year(1917).century, 20);
      assert.equal(new Year(1900).century, 19);
      assert.equal(new Year(2001).century, 21);
    });
  });

  describe('comparison', () => {
    it('compares by value, not by identity', () => {
      assert.ok(new Year(1969).equals(new Year(1969)));
      assert.ok(new Year(1969) !== new Year(1969)); // two distinct objects
    });

    it('orders years', () => {
      assert.ok(new Year(1900).isBefore(new Year(1901)));
      assert.ok(new Year(1901).isAfter(new Year(1900)));
    });
  });

  describe('range', () => {
    it('is inclusive at both ends', () => {
      const years = Year.range(new Year(1900), new Year(1903)).map((y) => y.value);
      assert.deepEqual(years, [1900, 1901, 1902, 1903]);
    });

    it('returns one year when start equals end', () => {
      assert.equal(Year.range(new Year(1969), new Year(1969)).length, 1);
    });

    it('refuses a backwards range', () => {
      assert.throws(() => Year.range(new Year(2000), new Year(1900)), ValidationError);
    });

    it('mvpRange runs from 1900 to this year', () => {
      const range = Year.mvpRange();
      assert.equal(range[0].value, 1900);
      assert.equal(range.at(-1).value, Year.latestAllowed());
    });
  });

  it('serialises to a plain number', () => {
    assert.deepEqual(JSON.parse(JSON.stringify({ y: new Year(1969) })), { y: 1969 });
  });

  it('keeps its value private and unwritable', () => {
    const year = new Year(1969);
    // #value is a genuinely private field: not enumerable, not reachable from outside.
    assert.deepEqual(Object.keys(year), []);
    assert.equal(year.value, 1969);
  });
});

describe('Country', () => {
  it('normalises the code to uppercase and trims the name', () => {
    const c = new Country({ id: 1, code: 'fra', name: '  France ', region: 'Europe' });
    assert.equal(c.code, 'FRA');
    assert.equal(c.name, 'France');
  });

  it('rejects an invalid code', () => {
    for (const code of ['FR', 'FRAN', '12A', '', null, 42]) {
      assert.throws(
        () => new Country({ id: 1, code, name: 'X', region: 'Europe' }),
        ValidationError,
        `should reject ${JSON.stringify(code)}`,
      );
    }
  });

  it('rejects a blank name', () => {
    assert.throws(
      () => new Country({ id: 1, code: 'FRA', name: '   ', region: 'Europe' }),
      ValidationError,
    );
  });

  it('defaults an unknown region rather than failing', () => {
    assert.equal(new Country({ id: 1, code: 'FRA', name: 'France', region: '' }).region, 'Unknown');
  });

  // Entities compare by IDENTITY — the defining difference from Year.
  it('treats two objects with the same id as the same country', () => {
    const a = new Country({ id: 7, code: 'FRA', name: 'France', region: 'Europe' });
    const b = new Country({ id: 7, code: 'FRA', name: 'French Republic', region: 'Europe' });
    assert.ok(a.equals(b));
  });

  it('is frozen after construction', () => {
    const c = new Country({ id: 1, code: 'FRA', name: 'France', region: 'Europe' });
    assert.throws(() => {
      c.name = 'Elsewhere';
    }, TypeError);
  });
});

describe('HistoricalEvent', () => {
  it('builds a valid event', () => {
    const e = anEvent({ year: 1969, title: 'Apollo 11' });
    assert.equal(e.year.value, 1969);
    assert.equal(e.title, 'Apollo 11');
  });

  // THE EDITORIAL RULE OF THE PROJECT, enforced in code.
  it('refuses to exist without a citable source', () => {
    for (const sourceUrl of ['', 'not-a-url', 'ftp://example.com/x', 'wikipedia.org/wiki/X', null]) {
      assert.throws(
        () => anEvent({ sourceUrl }),
        ValidationError,
        `should reject ${JSON.stringify(sourceUrl)}`,
      );
    }
  });

  it('requires a real Year, not a bare number', () => {
    assert.throws(
      () =>
        new HistoricalEvent({
          id: 1, countryId: 1, year: 1969, title: 'x', summary: 'y',
          category: 'politics', sourceUrl: 'https://example.org/a',
        }),
      ValidationError,
    );
  });

  it('rejects a title or summary over the limit', () => {
    assert.throws(() => anEvent({ title: 'x'.repeat(161) }), ValidationError);
    assert.throws(() => anEvent({ summary: 'x'.repeat(1001) }), ValidationError);
  });

  it('rejects an unknown category', () => {
    assert.throws(() => anEvent({ category: 'gossip' }), ValidationError);
  });

  it('rejects a malformed monthDay', () => {
    for (const monthDay of ['3-2', '0302', '03/02']) {
      assert.throws(() => anEvent({ monthDay }), ValidationError, `should reject ${monthDay}`);
    }
  });

  it('allows a year-only event', () => {
    assert.equal(anEvent({ monthDay: null }).monthDay, null);
  });

  describe('compareByDate', () => {
    it('orders by year first', () => {
      assert.ok(anEvent({ year: 1900 }).compareByDate(anEvent({ year: 1901 })) < 0);
    });

    it('orders by month and day within a year', () => {
      const march = anEvent({ year: 1969, monthDay: '03-02' });
      const april = anEvent({ year: 1969, monthDay: '04-28' });
      assert.ok(march.compareByDate(april) < 0);
    });

    it('puts dated events before undated ones', () => {
      const dated = anEvent({ year: 1969, monthDay: '03-02' });
      const undated = anEvent({ year: 1969, monthDay: null });
      assert.ok(dated.compareByDate(undated) < 0);
      assert.ok(undated.compareByDate(dated) > 0);
    });

    it('falls back to title for two undated events', () => {
      const a = anEvent({ year: 1969, monthDay: null, title: 'Alpha' });
      const b = anEvent({ year: 1969, monthDay: null, title: 'Beta' });
      assert.ok(a.compareByDate(b) < 0);
    });
  });

  it('serialises the year as a number for the API', () => {
    assert.equal(anEvent({ year: 1969 }).toJSON().year, 1969);
  });
});

describe('domain errors', () => {
  it('cannot construct the abstract base directly', () => {
    assert.throws(() => new DomainError('nope'), TypeError);
  });

  it('carries a machine code and an HTTP status', () => {
    const validation = new ValidationError('bad', 'year');
    assert.equal(validation.code, 'VALIDATION_ERROR');
    assert.equal(validation.httpStatus, 400);

    const notFound = new NotFoundError('Country', 'ZZZ');
    assert.equal(notFound.code, 'NOT_FOUND');
    assert.equal(notFound.httpStatus, 404);
    assert.match(notFound.message, /ZZZ/);
  });

  it('reports its own subclass name', () => {
    assert.equal(new ValidationError('x').name, 'ValidationError');
    assert.ok(new ValidationError('x') instanceof DomainError);
  });
});
