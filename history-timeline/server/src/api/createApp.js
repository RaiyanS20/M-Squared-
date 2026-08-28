import path from 'node:path';
import { fileURLToPath } from 'node:url';
import express from 'express';
import cors from 'cors';
import { TimelineController } from './controllers/TimelineController.js';
import { HealthController } from './controllers/HealthController.js';
import { errorHandler, notFoundHandler } from './errorHandler.js';

const here = path.dirname(fileURLToPath(import.meta.url));
const CLIENT_DIR = path.resolve(here, '../../../client');

/**
 * Builds the Express application from its dependencies. It does NOT listen on a
 * port — that is `main.js`'s job.
 *
 * Why the split? Because a test can call `createApp()` and drive it without
 * binding a socket, without a database, and without a free port. Separating
 * "build the app" from "run the app" is the single change that makes an HTTP
 * layer testable.
 *
 * This app also SERVES THE BROWSER APP as static files. With no build step and
 * no dev server, the browser loads everything from one origin — which means no
 * proxy configuration and no CORS in development. One less moving part to learn
 * while you are learning everything else.
 */
export function createApp({ service, db = null, corsOrigins = [], isProduction = false, serveClient = true }) {
  const app = express();

  // Trust the X-Forwarded-* headers set by an AWS load balancer, so req.ip is
  // the real client rather than the proxy.
  app.set('trust proxy', true);
  // Stop advertising "Express" to anyone scanning for known vulnerabilities.
  app.disable('x-powered-by');

  // Only needed when the browser app is served from a DIFFERENT origin than the
  // API (which is the production setup — see docs/10). Empty list = deny all,
  // never '*', which would tell every site on the internet it may read our
  // responses.
  if (corsOrigins.length > 0) {
    app.use(cors({ origin: corsOrigins, methods: ['GET'] }));
  }

  app.use(express.json({ limit: '100kb' }));

  app.use('/', new HealthController(db).routes());
  app.use('/api', new TimelineController(service).routes());

  if (serveClient) {
    app.use(express.static(CLIENT_DIR, { extensions: ['html'] }));
  }

  // ORDER MATTERS, and it is behaviour rather than style: 404 after every route,
  // and the error handler last of all. Put the error handler first and it never
  // fires; put notFoundHandler first and every request is a 404.
  app.use(notFoundHandler);
  app.use(errorHandler(isProduction));

  return app;
}
