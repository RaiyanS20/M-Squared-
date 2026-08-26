import 'dotenv/config';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Config } from '../config/index.js';
import { PostgresDatabase } from '../db/Database.js';
import { Migrator } from '../db/Migrator.js';

/**
 * `npm run db:migrate`
 *
 * In ESM there is no `__dirname`, so we derive it from `import.meta.url`. The
 * build copies the .sql files next to the compiled JS (see scripts/copyAssets.mjs)
 * so this same path works from `src/` under tsx and from `dist/` in production.
 */
const here = path.dirname(fileURLToPath(import.meta.url));
const migrationsDir = path.resolve(here, '../db/migrations');

async function main(): Promise<void> {
  const config = Config.fromEnv();
  const db = new PostgresDatabase({
    connectionString: config.databaseUrl,
    ssl: config.databaseSsl,
  });

  try {
    const applied = await new Migrator(db, migrationsDir).migrate();
    if (applied.length === 0) {
      console.log('[migrate] database is already up to date');
    } else {
      applied.forEach((f) => console.log(`[migrate] applied ${f}`));
    }
  } finally {
    await db.close();
  }
}

main().catch((error) => {
  console.error('[migrate] failed', error);
  process.exit(1);
});
