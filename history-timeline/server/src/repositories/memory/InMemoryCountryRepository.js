import { CountryRepository } from '../CountryRepository.js';

/**
 * A CountryRepository backed by plain arrays.
 *
 * OOP lesson: LISKOV SUBSTITUTION. This class can stand in for the Postgres one
 * anywhere, because it honours the same contract — same method names, same
 * return types, same ordering guarantees. If a test passes here and fails
 * against Postgres, one of the two broke the contract, and the shared suite in
 * `tests/repositoryContract.test.js` will say which.
 */
export class InMemoryCountryRepository extends CountryRepository {
  #data;

  constructor(dataset) {
    super();
    this.#data = dataset;
  }

  async findAll() {
    // Copy before sorting: `Array.prototype.sort` mutates in place, and this
    // array belongs to the dataset, not to us. Mutating shared data from a
    // "read" method is a classic and very confusing bug.
    return [...this.#data.countries].sort((a, b) => a.name.localeCompare(b.name));
  }

  async findAllWithEventsInYear(year) {
    const countryIds = new Set(
      this.#data.events.filter((e) => e.year.value === year).map((e) => e.countryId),
    );
    return (await this.findAll()).filter((c) => countryIds.has(c.id));
  }

  async findById(id) {
    return this.#data.countries.find((c) => c.id === id) ?? null;
  }

  async findByCode(code) {
    const upper = String(code).toUpperCase();
    return this.#data.countries.find((c) => c.code === upper) ?? null;
  }
}
