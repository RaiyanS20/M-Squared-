import type { NextFunction, Request, Response } from 'express';
import { DomainError } from '../domain/index.js';

/**
 * The single place HTTP status codes are chosen.
 *
 * OOP lesson: POLYMORPHISM paying rent. This function does not know about
 * `ValidationError` or `NotFoundError` specifically — it asks any `DomainError`
 * for its own `httpStatus` and `code`. Add a `ConflictError` tomorrow with
 * `httpStatus = 409` and this file does not change.
 *
 * Anything that is NOT a DomainError is, by definition, a bug we did not
 * anticipate: log it in full, and tell the client only "500", never the stack.
 * Leaking stack traces tells an attacker your file layout and library versions.
 */
export function errorHandler(isProduction: boolean) {
  return (error: unknown, _req: Request, res: Response, next: NextFunction): void => {
    if (res.headersSent) {
      next(error);
      return;
    }

    if (error instanceof DomainError) {
      res.status(error.httpStatus).json({
        error: { code: error.code, message: error.message },
      });
      return;
    }

    console.error('[api] unhandled error', error);
    res.status(500).json({
      error: {
        code: 'INTERNAL_ERROR',
        message: isProduction
          ? 'Something went wrong. Please try again.'
          : error instanceof Error
            ? error.message
            : String(error),
      },
    });
  };
}

/** 404 for unknown routes, so the client always gets JSON instead of HTML. */
export function notFoundHandler(req: Request, res: Response): void {
  res.status(404).json({
    error: { code: 'ROUTE_NOT_FOUND', message: `No route matches ${req.method} ${req.path}.` },
  });
}
