import { Router } from 'express';
import { ValidationError, Year } from '../../domain/index.js';
import { asyncRoute } from '../asyncRoute.js';

/**
 * Translates HTTP into service calls, and domain objects back into JSON.
 *
 * OOP lesson: keep controllers THIN. Every line here does one of exactly three
 * jobs: read input, call one service method, shape the response. The moment a
 * controller contains an `if` about history, that rule belongs in the service.
 *
 * The controller receives its service through the constructor rather than
 * importing a global — which is what lets `tests/api.test.js` boot the entire
 * HTTP stack against in-memory data.
 */
export class TimelineController {
  #service;

  constructor(service) {
    this.#service = service;
  }

  /** Builds the Express sub-router mounted at /api. */
  routes() {
    const router = Router();
    router.get('/timeline', asyncRoute(this.#getTimeline));
    router.get('/countries', asyncRoute(this.#listCountries));
    router.get('/years/:year/countries', asyncRoute(this.#listCountriesForYear));
    router.get('/years/:year/events', asyncRoute(this.#listEventsForYear));
    router.get('/years/:year/countries/:code/events', asyncRoute(this.#listEventsForYearAndCountry));
    return router;
  }

  /**
   * These are arrow-function fields, not methods, and it matters.
   *
   * A normal method loses its `this` when detached:
   *     router.get('/timeline', this.getTimeline)   // `this` is undefined inside
   * because Express calls the function without a receiver. Arrow fields capture
   * `this` at construction, so the reference is safe to pass around. (The older
   * fix is `.bind(this)` in the constructor — same effect, more noise.)
   *
   * This is the single most common "why is `this` undefined" bug in JavaScript.
   */
  #getTimeline = async (req, res) => {
    const start = this.#optionalYear(req.query.startYear, 'startYear');
    const end = this.#optionalYear(req.query.endYear, 'endYear');
    res.json(await this.#service.getTimelineScale(start, end));
  };

  #listCountries = async (_req, res) => {
    res.json({ countries: await this.#service.listCountries() });
  };

  #listCountriesForYear = async (req, res) => {
    const year = this.#requiredYear(req.params.year);
    res.json({ year, countries: await this.#service.listCountriesForYear(year) });
  };

  #listEventsForYear = async (req, res) => {
    const year = this.#requiredYear(req.params.year);
    res.json({ year, events: await this.#service.getEventsForYear(year) });
  };

  #listEventsForYearAndCountry = async (req, res) => {
    const year = this.#requiredYear(req.params.year);
    const code = String(req.params.code ?? '');
    // Validating the SHAPE ("is this three letters?") belongs here. Validating
    // the MEANING ("does this country exist?") belongs in the service.
    if (!/^[A-Za-z]{3}$/.test(code)) {
      throw new ValidationError(`"${code}" is not a 3-letter country code.`, 'code');
    }
    res.json(await this.#service.getEventsForYearAndCountry(year, code));
  };

  /** Path params are always strings. `Year.fromString` throws a 400-mapped error. */
  #requiredYear(raw) {
    return Year.fromString(raw).value;
  }

  #optionalYear(raw, field) {
    if (raw === undefined || raw === '') return undefined;
    // `?startYear=1&startYear=2` makes Express hand you an ARRAY. Without this
    // check that becomes a confusing failure deeper in the stack.
    if (typeof raw !== 'string') {
      throw new ValidationError(`${field} must be a single value.`, field);
    }
    return Year.fromString(raw).value;
  }
}
