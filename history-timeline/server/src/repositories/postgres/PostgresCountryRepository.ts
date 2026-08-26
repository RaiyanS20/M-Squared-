import { CountryRepository } from '../CountryRepository.js';
import type { Country } from '../../domain/index.js';
import type { QueryRunner } from '../../db/Database.js';
import { RowMapper, type CountryRow } from './rowMappers.js';

/**
 * The production CountryRepository.
 *
 * Every query uses PARAMETERISED placeholders ($1, $2). Never build SQL by
 * concatenating user input — that is SQL injection, and it is the single most
 * common way a project like this gets compromised. `pg` sends the SQL and the
 * values separately, so a value can never be parsed as SQL.
 */
export class PostgresCountryRepository extends CountryRepository {
  constructor(private readonly db: QueryRunner) {
    super();
  }

  async findAll(): Promise<Country[]> {
    const rows = await this.db.query<CountryRow>(
      'SELECT id, code, name, region FROM countries ORDER BY name ASC',
    );
    return rows.map(RowMapper.toCountry);
  }

  async findAllWithEventsInYear(year: number): Promise<Country[]> {
    const rows = await this.db.query<CountryRow>(
      `SELECT c.id, c.code, c.name, c.region
         FROM countries c
        WHERE EXISTS (
              SELECT 1 FROM historical_events e
               WHERE e.country_id = c.id AND e.year = $1
        )
        ORDER BY c.name ASC`,
      [year],
    );
    return rows.map(RowMapper.toCountry);
  }

  async findById(id: number): Promise<Country | null> {
    const rows = await this.db.query<CountryRow>(
      'SELECT id, code, name, region FROM countries WHERE id = $1',
      [id],
    );
    return rows[0] ? RowMapper.toCountry(rows[0]) : null;
  }

  async findByCode(code: string): Promise<Country | null> {
    const rows = await this.db.query<CountryRow>(
      'SELECT id, code, name, region FROM countries WHERE code = $1',
      [code.toUpperCase()],
    );
    return rows[0] ? RowMapper.toCountry(rows[0]) : null;
  }
}
