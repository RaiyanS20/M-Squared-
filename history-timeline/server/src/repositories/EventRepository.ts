import type { HistoricalEvent } from '../domain/index.js';

/** One row of the "how many events per year" histogram behind the timeline. */
export interface YearEventCount {
  year: number;
  eventCount: number;
}

/**
 * The contract every event data source must satisfy.
 *
 * OOP lesson: INTERFACE SEGREGATION. This repository knows about events and
 * nothing else. A fat `DataRepository` with 40 methods would force the
 * in-memory test double to implement things it never uses.
 */
export abstract class EventRepository {
  /** Events in one country in one year, oldest-dated first. */
  abstract findByYearAndCountry(year: number, countryId: number): Promise<HistoricalEvent[]>;

  /** Every event in a year, across all countries. */
  abstract findByYear(year: number): Promise<HistoricalEvent[]>;

  /** Event counts per year across an inclusive range — powers the timeline density bars. */
  abstract countByYearRange(startYear: number, endYear: number): Promise<YearEventCount[]>;
}
