import 'dotenv/config';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Config } from '../config/index.js';
import { PostgresDatabase } from '../db/Database.js';
import { Migrator } from '../db/Migrator.js';

/**
 * `npm run db:migrate`
 *
 * ES modules have no `__dirname`, so we derive it from `import.meta.url`. This
 * is one of the few places the module system is visibly different from the
 * older CommonJS (`require`) style you will still see in older tutorials.
 */
const here = path.dirname(fileURLToPath(import.meta.url));
const migrationsDir = path.resolve(here, '../db/migrations');

async function main() {
  const config = Config.fromEnv();
  const db = new PostgresDatabase({
    connectionString: config.databaseUrl,
    ssl: config.databaseSsl,
  });

  try {
    const applied = await new Migrator(db, migrationsDir).migrate();
    if (applied.length === 0) console.log('[migrate] database is already up to date');
    else applied.forEach((f) => console.log(`[migrate] applied ${f}`));
  } finally {
    // `finally` runs whether or not the migration threw, so the pool is always
    // closed and the process can actually exit.
    await db.close();
  }
}

main().catch((error) => {
  console.error('[migrate] failed', error);
  process.exit(1);
});
