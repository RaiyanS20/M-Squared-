import { NotFoundError, ValidationError, Year } from '../domain/index.js';

/**
 * The application's use cases — one method per user-visible action.
 *
 * OOP lesson: SINGLE RESPONSIBILITY, and why a service layer earns its keep.
 *
 * Controllers translate HTTP. Repositories translate SQL. This class is the only
 * place that knows what the PRODUCT does. Notice it never mentions `req`, `res`,
 * `pg` or `SELECT` — which is exactly why its tests run in milliseconds with no
 * server and no database.
 *
 * Notice the constructor too: it takes the two ABSTRACT repositories. This class
 * has no idea whether it is talking to PostgreSQL or to two arrays.
 */
export class TimelineService {
  #countries;
  #events;

  constructor(countryRepository, eventRepository) {
    this.#countries = countryRepository;
    this.#events = eventRepository;
  }

  /**
   * MVP FEATURE 1: the timeline scale, with an event count per year so the UI
   * can show which years are dense with history.
   *
   * Defaults to the MVP window (1900 → today). The domain still permits 100 AD,
   * so widening the timeline later is a change of DEFAULTS, not of architecture.
   */
  async getTimelineScale(startYear, endYear) {
    const start = new Year(startYear ?? Year.MVP_START);
    const end = new Year(endYear ?? Year.latestAllowed());

    if (start.isAfter(end)) {
      throw new ValidationError(
        `startYear (${start.value}) must not be after endYear (${end.value}).`,
        'startYear',
      );
    }

    // ONE query for the whole range, then fill the gaps in memory. The
    // alternative — one query per year — is 126 round trips to render one page.
    // That is the N+1 problem, and it is invisible on a laptop.
    const counts = await this.#events.countByYearRange(start.value, end.value);
    const byYear = new Map(counts.map((c) => [c.year, c.eventCount]));

    const ticks = Year.range(start, end).map((year) => ({
      year: year.value,
      eventCount: byYear.get(year.value) ?? 0,
    }));

    return { startYear: start.value, endYear: end.value, ticks };
  }

  /** MVP FEATURE 3: every country we hold, for the browse list. */
  async listCountries() {
    return this.#countries.findAll();
  }

  /**
   * MVP FEATURES 2 + 3: after picking a year, show only the countries that
   * actually have something recorded. A list of 195 countries where 190 lead to
   * an empty page is a bad experience — and note that this is a PRODUCT
   * decision, made here, not in the UI.
   */
  async listCountriesForYear(year) {
    const validYear = new Year(year);
    return this.#countries.findAllWithEventsInYear(validYear.value);
  }

  /**
   * MVP FEATURE 4: the payoff — what happened in this country, in this year.
   *
   * Takes an ISO alpha-3 code because codes are stable and readable in a URL:
   * /api/years/1969/countries/FRA/events beats .../countries/74/events.
   */
  async getEventsForYearAndCountry(year, countryCode) {
    const validYear = new Year(year);

    const country = await this.#countries.findByCode(countryCode);
    if (!country) {
      // The repository returns null for "not here". Deciding whether that is an
      // ERROR is a product judgement, so it belongs in the service.
      throw new NotFoundError('Country', countryCode);
    }

    const events = await this.#events.findByYearAndCountry(validYear.value, country.id);
    return { year: validYear.value, country, events };
  }

  /** Everything that happened in a year, worldwide. A "surprise me" view. */
  async getEventsForYear(year) {
    return this.#events.findByYear(new Year(year).value);
  }
}
