import { CountryRepository } from '../CountryRepository.js';
import { RowMapper } from './rowMappers.js';

/**
 * The production CountryRepository.
 *
 * Every query uses PARAMETERISED placeholders ($1, $2). Never build SQL by
 * concatenating input — that is SQL injection, and it is the most common way a
 * project like this gets compromised. `pg` sends the SQL and the values to the
 * server separately, so a value can never be parsed as SQL.
 */
export class PostgresCountryRepository extends CountryRepository {
  #db;

  constructor(db) {
    super();
    this.#db = db;
  }

  async findAll() {
    const rows = await this.#db.query(
      'SELECT id, code, name, region FROM countries ORDER BY name ASC',
    );
    return rows.map(RowMapper.toCountry);
  }

  async findAllWithEventsInYear(year) {
    // EXISTS stops at the first match per country and cannot produce
    // duplicates, so there is no DISTINCT to deduplicate afterwards.
    const rows = await this.#db.query(
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

  async findById(id) {
    const rows = await this.#db.query(
      'SELECT id, code, name, region FROM countries WHERE id = $1',
      [id],
    );
    return rows[0] ? RowMapper.toCountry(rows[0]) : null;
  }

  async findByCode(code) {
    const rows = await this.#db.query(
      'SELECT id, code, name, region FROM countries WHERE code = $1',
      [String(code).toUpperCase()],
    );
    return rows[0] ? RowMapper.toCountry(rows[0]) : null;
  }
}
