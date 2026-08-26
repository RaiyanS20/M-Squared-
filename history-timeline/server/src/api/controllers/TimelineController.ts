import { Router, type Request, type Response } from 'express';
import { ValidationError, Year } from '../../domain/index.js';
import type { TimelineService } from '../../services/index.js';
import { asyncRoute } from '../asyncRoute.js';

/**
 * Translates HTTP into service calls, and domain objects back into JSON.
 *
 * OOP lesson: keep controllers THIN. Every line here is one of exactly three
 * jobs: read input, call one service method, shape the response. The moment a
 * controller contains an `if` about history, that rule belongs in the service.
 *
 * The controller receives its service through the constructor rather than
 * importing a global. That is what makes `tests/integration/api.test.ts` able to
 * boot the whole HTTP stack against in-memory data.
 */
export class TimelineController {
  constructor(private readonly service: TimelineService) {}

  /** Builds the Express sub-router for everything under /api. */
  routes(): Router {
    const router = Router();
    router.get('/timeline', asyncRoute(this.getTimeline));
    router.get('/countries', asyncRoute(this.listCountries));
    router.get('/years/:year/countries', asyncRoute(this.listCountriesForYear));
    router.get('/years/:year/events', asyncRoute(this.listEventsForYear));
    router.get('/years/:year/countries/:code/events', asyncRoute(this.listEventsForYearAndCountry));
    return router;
  }

  /**
   * Arrow-function properties, not methods.
   * A plain method loses `this` when passed as `router.get(..., this.getTimeline)`.
   * Arrow properties capture `this` at construction, so the reference is safe to
   * detach. (The alternative is `.bind(this)` in the constructor — same effect,
   * more noise.)
   */
  private getTimeline = async (req: Request, res: Response): Promise<void> => {
    const start = this.optionalYear(req.query.startYear, 'startYear');
    const end = this.optionalYear(req.query.endYear, 'endYear');
    const scale = await this.service.getTimelineScale(start, end);
    res.json(scale);
  };

  private listCountries = async (_req: Request, res: Response): Promise<void> => {
    const countries = await this.service.listCountries();
    res.json({ countries });
  };

  private listCountriesForYear = async (req: Request, res: Response): Promise<void> => {
    const year = this.requiredYear(req.params.year);
    const countries = await this.service.listCountriesForYear(year);
    res.json({ year, countries });
  };

  private listEventsForYear = async (req: Request, res: Response): Promise<void> => {
    const year = this.requiredYear(req.params.year);
    const events = await this.service.getEventsForYear(year);
    res.json({ year, events });
  };

  private listEventsForYearAndCountry = async (req: Request, res: Response): Promise<void> => {
    const year = this.requiredYear(req.params.year);
    const code = String(req.params.code ?? '');
    if (!/^[A-Za-z]{3}$/.test(code)) {
      throw new ValidationError(`"${code}" is not a 3-letter country code.`, 'code');
    }
    const result = await this.service.getEventsForYearAndCountry(year, code);
    res.json(result);
  };

  /**
   * Path params arrive as strings (Express types them loosely enough to include
   * arrays, so we normalise). `Year.fromString` throws a ValidationError, which
   * the error handler maps to 400 — the controller never writes a status itself.
   */
  private requiredYear(raw: unknown): number {
    return Year.fromString(typeof raw === 'string' ? raw : String(raw ?? '')).value;
  }

  private optionalYear(raw: unknown, field: string): number | undefined {
    if (raw === undefined || raw === '') return undefined;
    if (typeof raw !== 'string') {
      throw new ValidationError(`${field} must be a single value.`, field);
    }
    return Year.fromString(raw).value;
  }
}
