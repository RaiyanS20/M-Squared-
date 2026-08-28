import { ValidationError } from './errors.js';

/**
 * A year on the timeline.
 *
 * OOP lesson: VALUE OBJECT + ENCAPSULATION + IMMUTABILITY.
 *
 * A plain `number` can be -4, or 99999, or NaN, or the string "1969". A `Year`
 * cannot: the only way to make one is this constructor, and the constructor
 * enforces the rules. Once a `Year` exists anywhere in the system, every layer
 * downstream can trust it without re-checking.
 *
 * That idea has a name — *parse, don't validate* — and it is the single most
 * useful habit in this course. Check once, at the boundary, then work with a
 * value that cannot be wrong.
 *
 * WHY THIS MATTERS MORE IN JAVASCRIPT: there is no compiler to catch
 * `getEvents("FRA", 1969)` with the arguments swapped. In a typed language the
 * type system does some of this work for you. Here, YOUR CONSTRUCTOR IS THE
 * TYPE SYSTEM. Every rule you do not write is a rule that does not exist.
 */
export class Year {
  /**
   * `#value` is a genuinely private field — not a convention, a language
   * feature. `year.#value` from outside this class is a syntax error, and it
   * does not show up in `Object.keys` or `JSON.stringify`. Combined with a
   * getter and no setter, the value is immutable after construction.
   */
  #value;

  /** The project's timeline goes back to 100 AD (the long-term goal). */
  static EARLIEST = 100;

  /** The MVP only renders 1900 onwards. */
  static MVP_START = 1900;

  constructor(value) {
    if (typeof value !== 'number' || !Number.isInteger(value)) {
      throw new ValidationError(`Year must be a whole number, received ${JSON.stringify(value)}.`, 'year');
    }
    if (value < Year.EARLIEST) {
      throw new ValidationError(`Year must be ${Year.EARLIEST} AD or later, received ${value}.`, 'year');
    }
    if (value > Year.latestAllowed()) {
      throw new ValidationError(
        `Year must not be in the future (latest allowed is ${Year.latestAllowed()}), received ${value}.`,
        'year',
      );
    }
    this.#value = value;
  }

  get value() {
    return this.#value;
  }

  /**
   * A named STATIC FACTORY. Useful at the HTTP boundary, where everything
   * arrives as a string. It documents intent better than an overloaded
   * constructor and keeps string-parsing out of the constructor's job.
   */
  static fromString(raw) {
    const trimmed = String(raw ?? '').trim();
    if (trimmed === '' || !/^-?\d+$/.test(trimmed)) {
      throw new ValidationError(`"${raw}" is not a valid year.`, 'year');
    }
    return new Year(Number.parseInt(trimmed, 10));
  }

  /** The current year in UTC, so tests do not drift across timezones. */
  static latestAllowed() {
    return new Date().getUTCFullYear();
  }

  /** Inclusive list of years, oldest first. Renders the timeline scale. */
  static range(start, end) {
    if (start.isAfter(end)) {
      throw new ValidationError(
        `Range start (${start.value}) must not be after range end (${end.value}).`,
        'range',
      );
    }
    const years = [];
    for (let y = start.value; y <= end.value; y += 1) years.push(new Year(y));
    return years;
  }

  /** The default MVP scale: 1900 → the current year. */
  static mvpRange() {
    return Year.range(new Year(Year.MVP_START), new Year(Year.latestAllowed()));
  }

  /** Which decade this belongs to: 1917 → 1910. Used for UI grouping. */
  get decade() {
    return Math.floor(this.#value / 10) * 10;
  }

  /** Which century, 1-based: 1917 → 20. */
  get century() {
    return Math.floor((this.#value - 1) / 100) + 1;
  }

  isAfter(other) {
    return this.#value > other.value;
  }

  isBefore(other) {
    return this.#value < other.value;
  }

  /**
   * Value objects compare by VALUE, not by reference.
   * `new Year(1969) === new Year(1969)` is false — two different objects.
   * `.equals()` is how you ask the question you actually meant.
   */
  equals(other) {
    return other instanceof Year && other.value === this.#value;
  }

  toString() {
    return String(this.#value);
  }

  /** Called automatically by JSON.stringify, so a Year serialises as a number. */
  toJSON() {
    return this.#value;
  }
}
