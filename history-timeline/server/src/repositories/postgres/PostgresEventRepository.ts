import { EventRepository, type YearEventCount } from '../EventRepository.js';
import type { HistoricalEvent } from '../../domain/index.js';
import type { QueryRunner } from '../../db/Database.js';
import { RowMapper, type EventRow } from './rowMappers.js';

const EVENT_COLUMNS = 'id, country_id, year, month_day, title, summary, category, source_url';

export class PostgresEventRepository extends EventRepository {
  constructor(private readonly db: QueryRunner) {
    super();
  }

  async findByYearAndCountry(year: number, countryId: number): Promise<HistoricalEvent[]> {
    const rows = await this.db.query<EventRow>(
      `SELECT ${EVENT_COLUMNS}
         FROM historical_events
        WHERE year = $1 AND country_id = $2
        ORDER BY month_day ASC NULLS LAST, title ASC`,
      [year, countryId],
    );
    return rows.map(RowMapper.toEvent);
  }

  async findByYear(year: number): Promise<HistoricalEvent[]> {
    const rows = await this.db.query<EventRow>(
      `SELECT ${EVENT_COLUMNS}
         FROM historical_events
        WHERE year = $1
        ORDER BY month_day ASC NULLS LAST, title ASC`,
      [year],
    );
    return rows.map(RowMapper.toEvent);
  }

  async countByYearRange(startYear: number, endYear: number): Promise<YearEventCount[]> {
    // COUNT(*) returns BIGINT, which `pg` hands back as a STRING to avoid
    // silently truncating values above 2^53. Parse it explicitly.
    const rows = await this.db.query<{ year: number; event_count: string }>(
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
