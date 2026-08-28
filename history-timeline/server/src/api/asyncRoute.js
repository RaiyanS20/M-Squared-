/**
 * Wraps an async handler so a rejected promise reaches the error handler.
 *
 * Express 5 forwards rejected promises automatically, but wrapping explicitly
 * keeps the intent visible — and keeps the code correct if it is ever
 * back-ported to Express 4, where an unhandled rejection makes the request hang
 * forever with no error at all.
 */
export function asyncRoute(handler) {
  return (req, res, next) => {
    handler(req, res).catch(next);
  };
}
