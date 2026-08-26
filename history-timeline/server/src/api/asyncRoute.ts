import type { NextFunction, Request, RequestHandler, Response } from 'express';

/**
 * Wraps an async handler so a rejected promise reaches Express's error handler.
 *
 * Express 5 forwards rejected promises automatically, but wrapping explicitly
 * keeps the intent obvious and keeps the code correct if it is ever back-ported
 * to Express 4, where an unhandled rejection silently hangs the request forever.
 */
export function asyncRoute(
  handler: (req: Request, res: Response) => Promise<void>,
): RequestHandler {
  return (req: Request, res: Response, next: NextFunction) => {
    handler(req, res).catch(next);
  };
}
