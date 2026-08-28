import { ApiClient } from './ApiClient.js';

/**
 * The one object the UI uses to reach the backend.
 *
 * Each method maps to exactly one MVP feature, and each is a single line,
 * because `ApiClient` already owns every HTTP concern.
 */
export class TimelineApi extends ApiClient {
  /** Feature 1: the years to render, with an event count for each. */
  getTimeline(startYear, endYear, signal) {
    return this.get(`/api/timeline${this.query({ startYear, endYear })}`, signal);
  }

  /** Feature 3: every country. */
  async getAllCountries(signal) {
    const { countries } = await this.get('/api/countries', signal);
    return countries;
  }

  /** Features 2 + 3: countries with something recorded in the chosen year. */
  getCountriesForYear(year, signal) {
    return this.get(`/api/years/${year}/countries`, signal);
  }

  /** Feature 4: the events themselves. */
  getEvents(year, countryCode, signal) {
    return this.get(`/api/years/${year}/countries/${countryCode}/events`, signal);
  }
}

/**
 * One shared instance for the app. The base URL is empty because the API serves
 * this page, so '/api/...' is same-origin.
 *
 * Tests construct their own `TimelineApi` against a stubbed `fetch` rather than
 * importing this, which is why the class and the instance are kept separate.
 */
export const timelineApi = new TimelineApi();
