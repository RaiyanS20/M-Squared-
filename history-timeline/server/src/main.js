import 'dotenv/config';
import { Config } from './config/index.js';
import { PostgresDatabase } from './db/Database.js';
import {
  PostgresCountryRepository,
  PostgresEventRepository,
} from './repositories/postgres/index.js';
import { TimelineService } from './services/index.js';
import { createApp } from './api/index.js';

/**
 * THE COMPOSITION ROOT.
 *
 * This is the one file allowed to know every concrete class in the system. It
 * wires the object graph once, at startup:
 *
 *   Config → PostgresDatabase → Repositories → TimelineService → Controllers
 *
 * Everything else receives what it needs through a constructor. That is why
 * TimelineService can be tested with arrays, and why swapping PostgreSQL for
 * something else would touch these fifteen lines and nothing more.
 */
async function bootstrap() {
  const config = Config.fromEnv();

  const db = new PostgresDatabase({
    connectionString: config.databaseUrl,
    ssl: config.databaseSsl,
  });

  const service = new TimelineService(
    new PostgresCountryRepository(db),
    new PostgresEventRepository(db),
  );

  const app = createApp({
    service,
    db,
    corsOrigins: config.corsOrigins,
    isProduction: config.isProduction,
  });

  const server = app.listen(config.port, () => {
    console.log(`[api] listening on http://localhost:${config.port} (${config.nodeEnv})`);
    console.log(`[api] open http://localhost:${config.port} in your browser`);
  });

  /**
   * GRACEFUL SHUTDOWN. When AWS ECS replaces this task it sends SIGTERM and then
   * waits. Without this handler the process is killed mid-request and users see
   * connection resets during every single deploy.
   */
  const shutdown = async (signal) => {
    console.log(`[api] ${signal} received, shutting down`);
    server.close(async () => {
      await db.close();
      process.exit(0);
    });
    // Do not wait forever for a stuck connection. `.unref()` stops this timer
    // itself from keeping the process alive.
    setTimeout(() => process.exit(1), 10_000).unref();
  };

  process.on('SIGTERM', () => void shutdown('SIGTERM'));
  process.on('SIGINT', () => void shutdown('SIGINT'));
}

bootstrap().catch((error) => {
  console.error('[api] failed to start', error);
  process.exit(1);
});
