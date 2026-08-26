import type { Country } from '../domain/index.js';

/**
 * The contract every country data source must satisfy.
 *
 * OOP lesson: ABSTRACTION + DEPENDENCY INVERSION (the "D" in SOLID).
 *
 * Note what is NOT here: no SQL, no `pg`, no HTTP. The service layer depends on
 * this abstract class, never on PostgreSQL. That inversion is what lets us run
 * the entire test suite against an in-memory implementation in milliseconds and
 * still ship Postgres to production — the code under test is literally the same.
 *
 * We use an `abstract class` rather than an `interface` because it survives
 * compilation: `instanceof` works, and it can hold shared helper logic later.
 */
export abstract class CountryRepository {
  /** Every country we know about, alphabetical by name. */
  abstract findAll(): Promise<Country[]>;

  /** Countries that have at least one recorded event in the given year. */
  abstract findAllWithEventsInYear(year: number): Promise<Country[]>;

  abstract findById(id: number): Promise<Country | null>;

  /** Lookup by ISO alpha-3 code, case-insensitive. */
  abstract findByCode(code: string): Promise<Country | null>;
}
