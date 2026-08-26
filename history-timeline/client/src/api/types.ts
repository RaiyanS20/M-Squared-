/**
 * The shapes the API returns.
 *
 * These mirror the `toJSON()` methods of the server's domain objects. They are
 * hand-written rather than imported from the server, and that is a deliberate
 * teaching point: the client and server are separately deployable programs, and
 * the JSON between them is a CONTRACT. Writing the contract down twice means a
 * server change that breaks it fails a client test instead of failing silently
 * in a user's browser.
 *
 * (In a larger project you would generate these from an OpenAPI schema — same
 * idea, automated. See docs/09.)
 */

export type EventCategory =
  | 'politics'
  | 'conflict'
  | 'science'
  | 'culture'
  | 'disaster'
  | 'economy'
  | 'society';

export interface Country {
  id: number;
  code: string;
  name: string;
  region: string;
}

export interface HistoricalEvent {
  id: number;
  countryId: number;
  year: number;
  title: string;
  summary: string;
  category: EventCategory;
  sourceUrl: string;
  monthDay: string | null;
}

export interface TimelineTick {
  year: number;
  eventCount: number;
}

export interface TimelineScale {
  startYear: number;
  endYear: number;
  ticks: TimelineTick[];
}

export interface CountriesForYear {
  year: number;
  countries: Country[];
}

export interface EventsForYearAndCountry {
  year: number;
  country: Country;
  events: HistoricalEvent[];
}
