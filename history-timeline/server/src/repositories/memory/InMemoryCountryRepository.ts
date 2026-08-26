import { CountryRepository } from '../CountryRepository.js';
import type { Country } from '../../domain/index.js';
import type { InMemoryDataset } from './InMemoryDataset.js';

/**
 * A CountryRepository backed by plain arrays.
 *
 * OOP lesson: LISKOV SUBSTITUTION. This class can stand in for the Postgres one
 * anywhere, because it honours the same contract — same method names, same
 * return types, same ordering guarantees. If a test passes here and fails
 * against Postgres, one of the two implementations broke the contract, and the
 * shared contract test suite (tests/unit/repositoryContract.test.ts) will say so.
 */
export class InMemoryCountryRepository extends CountryRepository {
  constructor(private readonly data: InMemoryDataset) {
    super();
  }

  async findAll(): Promise<Country[]> {
    return [...this.data.countries].sort((a, b) => a.name.localeCompare(b.name));
  }

  async findAllWithEventsInYear(year: number): Promise<Country[]> {
    const countryIds = new Set(
      this.data.events.filter((e) => e.year.value === year).map((e) => e.countryId),
    );
    return (await this.findAll()).filter((c) => countryIds.has(c.id));
  }

  async findById(id: number): Promise<Country | null> {
    return this.data.countries.find((c) => c.id === id) ?? null;
  }

  async findByCode(code: string): Promise<Country | null> {
    const upper = code.toUpperCase();
    return this.data.countries.find((c) => c.code === upper) ?? null;
  }
}
