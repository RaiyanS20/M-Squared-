import { EventRepository } from '../EventRepository.js';
import { RowMapper } from './rowMappers.js';

/**
 * A hard-coded constant, NOT user input. This is the only interpolation in any
 * SQL string in the project, and it is worth knowing the difference:
 * interpolating a value you wrote is fine; interpolating a value a user sent is
 * how databases get dumped.
 */
const EVENT_COLUMNS = 'id, country_id, year, month_day, title, summary, category, source_url';

export class PostgresEventRepository extends EventRepository {
  #db;

  constructor(db) {
    super();
    this.#db = db;
  }

  async findByYearAndCountry(year, countryId) {
    const rows = await this.#db.query(
      `SELECT ${EVENT_COLUMNS}
         FROM historical_events
        WHERE year = $1 AND country_id = $2
        ORDER BY month_day ASC NULLS LAST, title ASC`,
      [year, countryId],
    );
    return rows.map(RowMapper.toEvent);
  }

  async findByYear(year) {
    const rows = await this.#db.query(
      `SELECT ${EVENT_COLUMNS}
         FROM historical_events
        WHERE year = $1
        ORDER BY month_day ASC NULLS LAST, title ASC`,
      [year],
    );
    return rows.map(RowMapper.toEvent);
  }

  async countByYearRange(startYear, endYear) {
    // COUNT(*) returns BIGINT, which exceeds JavaScript's safe integer range, so
    // `pg` hands it back as a STRING rather than silently losing precision.
    // Parse it explicitly — a library choosing correctness over convenience.
    const rows = await this.#db.query(
      `SELECT year, COUNT(*) AS event_count
         FROM historical_events
        WHERE year BETWEEN $1 AND $2
        GROUP BY year
        ORDER BY year ASC`,
      [startYear, endYear],
    );
    return rows.map((r) => ({ year: r.year, eventCount: Number.parseInt(r.event_count, 10) }));
  }
}
