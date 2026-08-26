import { Router, type Request, type Response } from 'express';
import type { PostgresDatabase } from '../../db/Database.js';
import { asyncRoute } from '../asyncRoute.js';

/**
 * Liveness and readiness.
 *
 * The distinction matters once this runs on AWS:
 *   * /health  — "is the process up?" Cheap. The load balancer uses this.
 *   * /ready   — "can it serve traffic?" Checks the database too, so a task with
 *                a dead DB connection is pulled out of rotation instead of
 *                returning 500s to real users.
 */
export class HealthController {
  constructor(private readonly db: PostgresDatabase | null) {}

  routes(): Router {
    const router = Router();
    router.get('/health', (_req: Request, res: Response) => {
      res.json({ status: 'ok', uptimeSeconds: Math.round(process.uptime()) });
    });
    router.get(
      '/ready',
      asyncRoute(async (_req, res) => {
        const dbHealthy = this.db ? await this.db.isHealthy() : true;
        res.status(dbHealthy ? 200 : 503).json({
          status: dbHealthy ? 'ready' : 'degraded',
          database: this.db ? (dbHealthy ? 'up' : 'down') : 'not-configured',
        });
      }),
    );
    return router;
  }
}
