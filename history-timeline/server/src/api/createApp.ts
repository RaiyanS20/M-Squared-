import express, { type Express } from 'express';
import cors from 'cors';
import type { TimelineService } from '../services/index.js';
import type { HealthCheck } from '../db/Database.js';
import { TimelineController } from './controllers/TimelineController.js';
import { HealthController } from './controllers/HealthController.js';
import { errorHandler, notFoundHandler } from './errorHandler.js';

export interface AppDependencies {
  service: TimelineService;
  /** Anything that can answer "is the datastore usable?" — see HealthCheck. */
  db?: HealthCheck | null;
  corsOrigins: string[];
  isProduction: boolean;
}

/**
 * Builds the Express application from its dependencies. It does NOT listen on a
 * port — that is `main.ts`'s job.
 *
 * Why the split? Because a test can call `createApp()` and drive it with
 * supertest without binding a socket, without a database, and without a free
 * port. Separating "build the app" from "run the app" is the single change that
 * makes an HTTP layer testable.
 */
export function createApp(deps: AppDependencies): Express {
  const app = express();

  // Trust the X-Forwarded-* headers set by an AWS Application Load Balancer, so
  // req.ip and req.protocol reflect the real client rather than the proxy.
  app.set('trust proxy', true);
  app.disable('x-powered-by');

  app.use(
    cors({
      origin: deps.corsOrigins.length > 0 ? deps.corsOrigins : false,
      methods: ['GET'],
    }),
  );
  app.use(express.json({ limit: '100kb' }));

  app.use('/', new HealthController(deps.db ?? null).routes());
  app.use('/api', new TimelineController(deps.service).routes());

  // Order matters: 404 first, then the error handler, and both AFTER every route.
  app.use(notFoundHandler);
  app.use(errorHandler(deps.isProduction));

  return app;
}
