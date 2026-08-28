import 'dotenv/config';
import { Config } from '../config/index.js';
import { PostgresDatabase } from '../db/Database.js';
import { SEED_COUNTRIES, SEED_EVENTS } from '../db/seedData.js';
import { Country, HistoricalEvent, Year } from '../domain/index.js';

/**
 * `npm run db:seed`
 *
 * Two properties make this safe to run repeatedly:
 *
 *   1. IDEMPOTENT — `ON CONFLICT ... DO UPDATE` means running it twice leaves
 *      the same rows, not duplicates. Correct a typo in seedData.js, re-run, and
 *      the fix lands. Any script you might run twice should be safe to run
 *      twice; it is one of the cheapest reliability habits there is.
 *
 *   2. VALIDATED FIRST — every row goes through the domain constructors BEFORE
 *      anything is written. A bad source URL or a 2,000-character summary fails
 *      here, loudly, instead of quietly entering the dataset.
 *
 * The write is one transaction, so a failure halfway leaves the database exactly
 * as it was.
 */
async function main() {
  const config = Config.fromEnv();
  const db = new PostgresDatabase({
    connectionString: config.databaseUrl,
    ssl: config.databaseSsl,
  });

  try {
    // --- Validate ---------------------------------------------------------
    // id 0 is a placeholder: ids are assigned by Postgres, and validation does
    // not depend on them.
    SEED_COUNTRIES.forEach((c) => new Country({ id: 0, ...c }));

    const countryCodes = new Set(SEED_COUNTRIES.map((c) => c.code));
    for (const e of SEED_EVENTS) {
      if (!countryCodes.has(e.countryCode)) {
        throw new Error(`Event "${e.title}" references unknown country ${e.countryCode}`);
      }
      new HistoricalEvent({ id: 0, countryId: 0, ...e, year: new Year(e.year) });
    }
    console.log(
      `[seed] validated ${SEED_COUNTRIES.length} countries and ${SEED_EVENTS.length} events`,
    );

    // --- Write ------------------------------------------------------------
    await db.transaction(async (runner) => {
      const idByCode = new Map();

      for (const c of SEED_COUNTRIES) {
        const rows = await runner.query(
          `INSERT INTO countries (code, name, region)
                VALUES ($1, $2, $3)
           ON CONFLICT (code) DO UPDATE SET name = EXCLUDED.name, region = EXCLUDED.region
             RETURNING id`,
          [c.code, c.name, c.region],
        );
        idByCode.set(c.code, rows[0].id);
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
          [idByCode.get(e.countryCode), e.year, e.monthDay, e.title, e.summary, e.category, e.sourceUrl],
        );
      }
    });

    const rows = await db.query('SELECT COUNT(*) AS count FROM historical_events');
    console.log(`[seed] done — ${rows[0]?.count ?? 0} events in the database`);
  } finally {
    await db.close();
  }
}

main().catch((error) => {
  console.error('[seed] failed', error);
  process.exit(1);
});
