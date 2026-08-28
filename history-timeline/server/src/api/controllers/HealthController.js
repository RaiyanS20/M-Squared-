import { Router } from 'express';
import { asyncRoute } from '../asyncRoute.js';

/**
 * Liveness and readiness.
 *
 * The distinction matters once this runs on AWS:
 *   /health — "is the process alive?" Cheap. Answers "should this be restarted?"
 *   /ready  — "can it serve traffic?" Checks the database too, so a task with a
 *             dead connection is pulled out of the load balancer instead of
 *             returning 500s to real users.
 *
 * `db` is anything with an `isHealthy()` method — not specifically a
 * PostgresDatabase. The API layer therefore names no storage type at all, and a
 * test can pass `null` or a two-line fake.
 */
export class HealthController {
  #db;

  constructor(db) {
    this.#db = db ?? null;
  }

  routes() {
    const router = Router();

    router.get('/health', (_req, res) => {
      res.json({ status: 'ok', uptimeSeconds: Math.round(process.uptime()) });
    });

    router.get(
      '/ready',
      asyncRoute(async (_req, res) => {
        const dbHealthy = this.#db ? await this.#db.isHealthy() : true;
        res.status(dbHealthy ? 200 : 503).json({
          status: dbHealthy ? 'ready' : 'degraded',
          database: this.#db ? (dbHealthy ? 'up' : 'down') : 'not-configured',
        });
      }),
    );

    return router;
  }
}
