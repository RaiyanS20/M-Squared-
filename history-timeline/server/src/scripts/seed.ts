import 'dotenv/config';
import { Config } from '../config/index.js';
import { PostgresDatabase, type QueryRunner } from '../db/Database.js';
import { SEED_COUNTRIES, SEED_EVENTS } from '../db/seedData.js';
import { Country, HistoricalEvent, Year } from '../domain/index.js';

/**
 * `npm run db:seed`
 *
 * Two properties make this script safe to run repeatedly:
 *
 *   1. IDEMPOTENT. `ON CONFLICT ... DO UPDATE` means running it twice leaves the
 *      same rows, not duplicates. You can re-run it after editing seedData.ts and
 *      the corrections land.
 *   2. VALIDATED. Every row is passed through the domain constructors BEFORE it
 *      reaches the database. A typo in a source URL or a 2000-character summary
 *      fails here, loudly, instead of silently entering the dataset.
 *
 * The whole thing runs in one transaction, so a failure halfway leaves the
 * database exactly as it was.
 */
async function main(): Promise<void> {
  const config = Config.fromEnv();
  const db = new PostgresDatabase({
    connectionString: config.databaseUrl,
    ssl: config.databaseSsl,
  });

  try {
    // --- Validate first, write second. -------------------------------------
    // Constructing the domain objects proves the data is sound. Note we pass a
    // placeholder id of 0: ids are assigned by Postgres, and validation does not
    // depend on them.
    SEED_COUNTRIES.forEach((c) => new Country({ id: 0, ...c }));

    const countryCodes = new Set(SEED_COUNTRIES.map((c) => c.code));
    for (const e of SEED_EVENTS) {
      if (!countryCodes.has(e.countryCode)) {
        throw new Error(`Event "${e.title}" references unknown country ${e.countryCode}`);
      }
      new HistoricalEvent({
        id: 0,
        countryId: 0,
        year: new Year(e.year),
        title: e.title,
        summary: e.summary,
        category: e.category,
        sourceUrl: e.sourceUrl,
        monthDay: e.monthDay,
      });
    }
    console.log(
      `[seed] validated ${SEED_COUNTRIES.length} countries and ${SEED_EVENTS.length} events`,
    );

    // --- Write. ------------------------------------------------------------
    await db.transaction(async (runner: QueryRunner) => {
      const idByCode = new Map<string, number>();

      for (const c of SEED_COUNTRIES) {
        const rows = await runner.query<{ id: number }>(
          `INSERT INTO countries (code, name, region)
                VALUES ($1, $2, $3)
           ON CONFLICT (code) DO UPDATE SET name = EXCLUDED.name, region = EXCLUDED.region
             RETURNING id`,
          [c.code, c.name, c.region],
        );
        idByCode.set(c.code, rows[0]!.id);
      }

      for (const e of SEED_EVENTS) {
        await runner.query(
          `INSERT INTO historical_events
                  (country_id, year, month_day, title, summary, category, source_url)
                VALUES ($1, $2, $3, $4, $5, $6, $7)
           ON CONFLICT (country_id, year, title) DO UPDATE
                   SET summary    = EXCLUDED.summary,
                       category   = EXCLUDED.category,
                       source_url = EXCLUDED.source_url,
                       month_day  = EXCLUDED.month_day`,
          [
            idByCode.get(e.countryCode),
            e.year,
            e.monthDay,
            e.title,
            e.summary,
            e.category,
            e.sourceUrl,
          ],
        );
      }
    });

    const rows = await db.query<{ count: string }>(
      'SELECT COUNT(*) AS count FROM historical_events',
    );
    console.log(`[seed] done — ${rows[0]?.count ?? '0'} events in the database`);
  } finally {
    await db.close();
  }
}

main().catch((error) => {
  console.error('[seed] failed', error);
  process.exit(1);
});
