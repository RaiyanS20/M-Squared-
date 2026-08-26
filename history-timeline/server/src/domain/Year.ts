import { ValidationError } from './errors.js';

/**
 * A year on the timeline.
 *
 * OOP lesson: VALUE OBJECT + ENCAPSULATION + IMMUTABILITY.
 *
 * A plain `number` can be -4 or 99_999 or NaN. A `Year` cannot: the only way to
 * build one is through the constructor, and the constructor enforces the rules.
 * Once a `Year` exists anywhere in the system, every layer downstream can trust
 * it without re-checking. This is "parse, don't validate" — push the check to
 * the boundary once, then work with a type that makes bad states unrepresentable.
 *
 * The field is `readonly`, so a Year is immutable: operations like `next()`
 * return a NEW Year rather than mutating this one. Immutable objects are safe to
 * share, cache and compare.
 */
export class Year {
  /** The project's timeline starts at 100 AD (the long-term goal). */
  static readonly EARLIEST = 100;

  /** The MVP timeline only renders 1900 onwards. */
  static readonly MVP_START = 1900;

  readonly value: number;

  constructor(value: number) {
    if (!Number.isInteger(value)) {
      throw new ValidationError(`Year must be a whole number, received "${value}".`, 'year');
    }
    if (value < Year.EARLIEST) {
      throw new ValidationError(
        `Year must be ${Year.EARLIEST} AD or later, received ${value}.`,
        'year',
      );
    }
    if (value > Year.latestAllowed()) {
      throw new ValidationError(
        `Year must not be in the future (latest allowed is ${Year.latestAllowed()}), received ${value}.`,
        'year',
      );
    }
    this.value = value;
    Object.freeze(this);
  }

  /**
   * Static factory. Useful at the HTTP boundary where everything is a string.
   * OOP lesson: a named factory documents intent better than an overloaded
   * constructor, and keeps the string-parsing concern out of the constructor.
   */
  static fromString(raw: string): Year {
    const trimmed = raw.trim();
    if (trimmed === '' || !/^-?\d+$/.test(trimmed)) {
      throw new ValidationError(`"${raw}" is not a valid year.`, 'year');
    }
    return new Year(Number.parseInt(trimmed, 10));
  }

  /** The current year, in UTC, so tests do not drift across timezones. */
  static latestAllowed(): number {
    return new Date().getUTCFullYear();
  }

  /** Inclusive list of years, oldest first. Used to render the timeline scale. */
  static range(start: Year, end: Year): Year[] {
    if (start.isAfter(end)) {
      throw new ValidationError(
        `Range start (${start.value}) must not be after range end (${end.value}).`,
        'range',
      );
    }
    const years: Year[] = [];
    for (let y = start.value; y <= end.value; y += 1) {
      years.push(new Year(y));
    }
    return years;
  }

  /** The default MVP scale: 1900 → the current year. */
  static mvpRange(): Year[] {
    return Year.range(new Year(Year.MVP_START), new Year(Year.latestAllowed()));
  }

  /** Which decade this year belongs to, e.g. 1917 → 1910. Used for UI grouping. */
  get decade(): number {
    return Math.floor(this.value / 10) * 10;
  }

  /** Which century this year belongs to, 1-based, e.g. 1917 → 20. */
  get century(): number {
    return Math.floor((this.value - 1) / 100) + 1;
  }

  isAfter(other: Year): boolean {
    return this.value > other.value;
  }

  isBefore(other: Year): boolean {
    return this.value < other.value;
  }

  /**
   * Value objects compare by VALUE, not by reference.
   * `new Year(1969) === new Year(1969)` is false; `.equals()` is how you ask.
   */
  equals(other: Year): boolean {
    return other instanceof Year && other.value === this.value;
  }

  toString(): string {
    return String(this.value);
  }

  /** Called automatically by `JSON.stringify`, so a Year serialises as a number. */
  toJSON(): number {
    return this.value;
  }
}
