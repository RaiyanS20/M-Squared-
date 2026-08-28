import { DomainError } from '../domain/index.js';

/**
 * The single place HTTP status codes are chosen.
 *
 * OOP lesson: POLYMORPHISM paying rent. This function does not know about
 * `ValidationError` or `NotFoundError` specifically — it asks any `DomainError`
 * for its own `httpStatus` and `code`. Add a `ConflictError` with
 * `httpStatus = 409` tomorrow and this file does not change. That is the
 * open/closed principle in production code rather than in a textbook.
 *
 * Anything that is NOT a DomainError is, by definition, a bug we did not
 * anticipate: log it in full, and tell the client only "500" — never the stack.
 * A stack trace tells an attacker your file layout and library versions.
 */
export function errorHandler(isProduction) {
  // Express identifies error middleware by ARITY: a function of four arguments.
  // Write three and it is registered as ordinary middleware that never sees an
  // error. This is a genuinely common and baffling bug.
  return (error, _req, res, next) => {
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

/** 404 for unknown API routes, so the client always gets JSON, never HTML. */
export function notFoundHandler(req, res) {
  res.status(404).json({
    error: { code: 'ROUTE_NOT_FOUND', message: `No route matches ${req.method} ${req.path}.` },
  });
}
