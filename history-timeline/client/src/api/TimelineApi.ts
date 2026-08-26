import { ApiClient } from './ApiClient.js';
import type {
  CountriesForYear,
  Country,
  EventsForYearAndCountry,
  TimelineScale,
} from './types.js';

/**
 * The one object the UI uses to reach the backend.
 *
 * Each method maps to exactly one MVP feature, and each is a single typed line
 * because `ApiClient` already owns every HTTP concern. No component ever calls
 * `fetch` directly — so changing authentication, retries or the base URL is a
 * change to one file.
 */
export class TimelineApi extends ApiClient {
  constructor(baseUrl: string = import.meta.env.VITE_API_BASE_URL ?? '') {
    super(baseUrl);
  }

  /** Feature 1: the years to render, with an event count for each. */
  getTimeline(startYear?: number, endYear?: number, signal?: AbortSignal): Promise<TimelineScale> {
    return this.get<TimelineScale>(`/api/timeline${this.query({ startYear, endYear })}`, signal);
  }

  /** Feature 3: every country. */
  async getAllCountries(signal?: AbortSignal): Promise<Country[]> {
    const { countries } = await this.get<{ countries: Country[] }>('/api/countries', signal);
    return countries;
  }

  /** Features 2 + 3: countries that have something recorded in the chosen year. */
  getCountriesForYear(year: number, signal?: AbortSignal): Promise<CountriesForYear> {
    return this.get<CountriesForYear>(`/api/years/${year}/countries`, signal);
  }

  /** Feature 4: the events themselves. */
  getEvents(year: number, countryCode: string, signal?: AbortSignal): Promise<EventsForYearAndCountry> {
    return this.get<EventsForYearAndCountry>(
      `/api/years/${year}/countries/${countryCode}/events`,
      signal,
    );
  }
}

/**
 * One shared instance for the app. Tests construct their own `TimelineApi`
 * against a stubbed `fetch` instead of importing this, which is why the class
 * and the instance are kept separate.
 */
export const timelineApi = new TimelineApi();
