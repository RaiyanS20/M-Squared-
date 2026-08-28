import { Country, HistoricalEvent, Year } from '../../domain/index.js';

/**
 * Translates database rows into domain objects.
 *
 * OOP lesson: ANTI-CORRUPTION LAYER.
 *
 * The database speaks snake_case and knows nothing about `Year`. The domain
 * speaks camelCase and refuses invalid data. This module is the ONLY place the
 * two vocabularies meet — so renaming a column touches one file, not forty.
 *
 * It is also where a real, boring bug gets fixed once: `CHAR(3)` is blank-padded
 * by PostgreSQL, so 'FRA' can come back as 'FRA '. `.trim()` here means no other
 * line of code ever has to know that.
 */
export class RowMapper {
  static toCountry(row) {
    return new Country({
      id: row.id,
      code: row.code.trim(),
      name: row.name,
      region: row.region,
    });
  }

  static toEvent(row) {
    return new HistoricalEvent({
      id: row.id,
      countryId: row.country_id,
      year: new Year(row.year),
      title: row.title,
      summary: row.summary,
      category: row.category,
      sourceUrl: row.source_url,
      monthDay: row.month_day?.trim() || null,
    });
  }
}
