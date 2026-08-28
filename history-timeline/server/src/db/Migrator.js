import { readdir, readFile } from 'node:fs/promises';
import path from 'node:path';

/**
 * Applies versioned `.sql` files in filename order, exactly once each.
 *
 * Why hand-roll this instead of installing a library? Because the whole
 * mechanism is about thirty lines — a ledger table, a sorted file list, a
 * transaction per file — and understanding it makes every migration tool you
 * meet later (Flyway, Prisma Migrate, Knex) obvious rather than magical.
 *
 * The rules that keep a team out of trouble:
 *   1. Never edit a migration that has been applied anywhere but your laptop.
 *      Write a new one instead.
 *   2. Each file runs in a transaction, so a failure leaves no half-state.
 *   3. Zero-pad filenames (001_, 002_) so lexical order is real order.
 */
export class Migrator {
  #db;
  #migrationsDir;

  constructor(db, migrationsDir) {
    this.#db = db;
    this.#migrationsDir = migrationsDir;
  }

  async #ensureLedger() {
    await this.#db.query(`
      CREATE TABLE IF NOT EXISTS schema_migrations (
        filename   TEXT PRIMARY KEY,
        applied_at TIMESTAMPTZ NOT NULL DEFAULT now()
      )
    `);
  }

  async #appliedFilenames() {
    const rows = await this.#db.query('SELECT filename FROM schema_migrations');
    return new Set(rows.map((r) => r.filename));
  }

  /** Returns the filenames applied by this run (empty if already up to date). */
  async migrate() {
    await this.#ensureLedger();
    const applied = await this.#appliedFilenames();

    const all = (await readdir(this.#migrationsDir)).filter((f) => f.endsWith('.sql')).sort();
    const pending = all.filter((f) => !applied.has(f));
    const executed = [];

    for (const filename of pending) {
      const sql = await readFile(path.join(this.#migrationsDir, filename), 'utf8');
      await this.#db.transaction(async (runner) => {
        await runner.query(sql);
        await runner.query('INSERT INTO schema_migrations (filename) VALUES ($1)', [filename]);
      });
      executed.push(filename);
    }

    return executed;
  }
}
