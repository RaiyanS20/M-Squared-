import {
  Country,
  HistoricalEvent,
  NotFoundError,
  ValidationError,
  Year,
} from '../domain/index.js';
import type { CountryRepository, EventRepository } from '../repositories/index.js';

/** One tick on the rendered timeline. */
export interface TimelineTick {
  year: number;
  eventCount: number;
}

export interface TimelineScale {
  startYear: number;
  endYear: number;
  ticks: TimelineTick[];
}

export interface YearInCountry {
  year: number;
  country: Country;
  events: HistoricalEvent[];
}

/**
 * The application's use cases, one method per user-visible action.
 *
 * OOP lesson: SINGLE RESPONSIBILITY, and the value of a service layer.
 * Controllers translate HTTP; repositories translate SQL; this class is the only
 * place that knows what the *product* does. Notice it never mentions `req`,
 * `res`, `pg` or `SELECT` — which is exactly why every test below runs in
 * milliseconds with no server and no database.
 *
 * Notice also the constructor: it accepts the two ABSTRACT repositories. This
 * class has no idea whether it is talking to Postgres or to two arrays.
 */
export class TimelineService {
  constructor(
    private readonly countries: CountryRepository,
    private readonly events: EventRepository,
  ) {}

  /**
   * Feature 1: the timeline scale, with an event count per year so the UI can
   * show which years are dense with history.
   *
   * Defaults to the MVP window (1900 → today). The domain still permits 100 AD,
   * so widening the MVP later is a change of DEFAULTS, not of architecture.
   */
  async getTimelineScale(startYear?: number, endYear?: number): Promise<TimelineScale> {
    const start = new Year(startYear ?? Year.MVP_START);
    const end = new Year(endYear ?? Year.latestAllowed());

    if (start.isAfter(end)) {
      throw new ValidationError(
        `startYear (${start.value}) must not be after endYear (${end.value}).`,
        'startYear',
      );
    }

    // One query for the whole range, then fill the gaps in memory. The
    // alternative — one query per year — is 126 round trips to render one page.
    const counts = await this.events.countByYearRange(start.value, end.value);
    const byYear = new Map(counts.map((c) => [c.year, c.eventCount]));

    const ticks = Year.range(start, end).map((year) => ({
      year: year.value,
      eventCount: byYear.get(year.value) ?? 0,
    }));

    return { startYear: start.value, endYear: end.value, ticks };
  }

  /** Feature 3: every country we hold, for the browse list. */
  async listCountries(): Promise<Country[]> {
    return this.countries.findAll();
  }

  /**
   * Feature 2 + 3 combined: after picking a year, show only the countries that
   * actually have something recorded — a list of 195 countries where 190 lead to
   * an empty page is a bad experience.
   */
  async listCountriesForYear(year: number): Promise<Country[]> {
    const validYear = new Year(year);
    return this.countries.findAllWithEventsInYear(validYear.value);
  }

  /**
   * Feature 4: the payoff — what happened in this country, in this year.
   * Accepts an ISO alpha-3 code because codes are stable and readable in a URL
   * (`/api/timeline/1969/FRA` beats `/api/timeline/1969/74`).
   */
  async getEventsForYearAndCountry(year: number, countryCode: string): Promise<YearInCountry> {
    const validYear = new Year(year);

    const country = await this.countries.findByCode(countryCode);
    if (!country) {
      throw new NotFoundError('Country', countryCode);
    }

    const events = await this.events.findByYearAndCountry(validYear.value, country.id);
    return { year: validYear.value, country, events };
  }

  /** Everything that happened in a year, worldwide. Useful for a "surprise me" view. */
  async getEventsForYear(year: number): Promise<HistoricalEvent[]> {
    return this.events.findByYear(new Year(year).value);
  }
}
