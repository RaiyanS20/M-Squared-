import { Country, HistoricalEvent, Year, EventCategory } from '../../domain/index.js';

/** Exactly the columns `countries` has, in snake_case as Postgres returns them. */
export interface CountryRow {
  id: number;
  code: string;
  name: string;
  region: string;
}

export interface EventRow {
  id: number;
  country_id: number;
  year: number;
  month_day: string | null;
  title: string;
  summary: string;
  category: string;
  source_url: string;
}

/**
 * Translates database rows into domain objects.
 *
 * OOP lesson: ANTI-CORRUPTION LAYER. The database speaks snake_case and knows
 * nothing about `Year` or `EventCategory`. The domain speaks camelCase and
 * refuses invalid data. This module is the only place the two vocabularies
 * meet — so a column rename touches one file, not forty.
 */
export class RowMapper {
  static toCountry(row: CountryRow): Country {
    return new Country({
      id: row.id,
      code: row.code.trim(), // CHAR(3) is blank-padded by Postgres.
      name: row.name,
      region: row.region,
    });
  }

  static toEvent(row: EventRow): HistoricalEvent {
    return new HistoricalEvent({
      id: row.id,
      countryId: row.country_id,
      year: new Year(row.year),
      title: row.title,
      summary: row.summary,
      category: row.category as EventCategory,
      sourceUrl: row.source_url,
      monthDay: row.month_day?.trim() || null,
    });
  }
}
