/**
 * The contract every event data source must satisfy.
 *
 * OOP lesson: INTERFACE SEGREGATION. This class knows about events and nothing
 * else. A fat `DataRepository` with forty methods would force the in-memory test
 * double to implement things it never uses.
 *
 * `countByYearRange` returns objects shaped `{ year, eventCount }`.
 */
export class EventRepository {
  constructor() {
    if (new.target === EventRepository) {
      throw new TypeError('EventRepository is abstract; extend it.');
    }
  }

  /** Events in one country in one year, oldest-dated first. */
  async findByYearAndCountry(_year, _countryId) {
    throw new Error(`${this.constructor.name} must implement findByYearAndCountry()`);
  }

  /** Every event in a year, across all countries. */
  async findByYear(_year) {
    throw new Error(`${this.constructor.name} must implement findByYear()`);
  }

  /** Event counts per year across an inclusive range — the timeline bars. */
  async countByYearRange(_startYear, _endYear) {
    throw new Error(`${this.constructor.name} must implement countByYearRange()`);
  }
}
