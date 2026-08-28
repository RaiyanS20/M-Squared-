/**
 * The contract every country data source must satisfy.
 *
 * OOP lesson: ABSTRACTION + DEPENDENCY INVERSION (the "D" in SOLID).
 *
 * Note what is NOT here: no SQL, no `pg`, no HTTP. The service layer depends on
 * this class, never on PostgreSQL. That inversion is why the whole test suite
 * can run against in-memory arrays in milliseconds while production runs on a
 * real database — the code under test is literally the same code.
 *
 * JavaScript has no `abstract` keyword and no `interface`. The pattern below is
 * how you express one anyway:
 *   * the constructor refuses to build the base class directly, and
 *   * each method throws unless a subclass overrides it.
 *
 * A typed language would catch a missing method at compile time. Here you find
 * out the first time it is called — which is exactly why the contract test in
 * `tests/repositoryContract.test.js` exists: it calls every method on every
 * implementation, so "the first time it is called" is in CI, not in production.
 */
export class CountryRepository {
  constructor() {
    if (new.target === CountryRepository) {
      throw new TypeError('CountryRepository is abstract; extend it.');
    }
  }

  /** Every country we know about, sorted by name. */
  async findAll() {
    throw new Error(`${this.constructor.name} must implement findAll()`);
  }

  /** Countries with at least one recorded event in the given year. */
  async findAllWithEventsInYear(_year) {
    throw new Error(`${this.constructor.name} must implement findAllWithEventsInYear()`);
  }

  async findById(_id) {
    throw new Error(`${this.constructor.name} must implement findById()`);
  }

  /** Lookup by ISO alpha-3 code, case-insensitive. */
  async findByCode(_code) {
    throw new Error(`${this.constructor.name} must implement findByCode()`);
  }
}
