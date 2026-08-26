import { readdir, readFile } from 'node:fs/promises';
import path from 'node:path';
import type { PostgresDatabase } from './Database.js';

/**
 * Applies versioned `.sql` files in filename order, exactly once each.
 *
 * Why hand-roll this instead of installing a migration library? Because the
 * whole mechanism is about thirty lines, and understanding it — a ledger table,
 * a sorted file list, a transaction per file — makes every migration tool you
 * meet later (Flyway, Prisma Migrate, Knex) obvious rather than magical.
 *
 * The rules that make migrations safe in a team:
 *   1. Never edit a migration that has been applied anywhere but your laptop.
 *      Write a new one instead.
 *   2. Each file runs inside a transaction, so a failure leaves no half-state.
 *   3. Filenames are zero-padded (`001_`, `002_`) so lexical order == real order.
 */
export class Migrator {
  constructor(
    private readonly db: PostgresDatabase,
    private readonly migrationsDir: string,
  ) {}

  private async ensureLedger(): Promise<void> {
    await this.db.query(`
      CREATE TABLE IF NOT EXISTS schema_migrations (
        filename   TEXT PRIMARY KEY,
        applied_at TIMESTAMPTZ NOT NULL DEFAULT now()
      )
    `);
  }

  private async appliedFilenames(): Promise<Set<string>> {
    const rows = await this.db.query<{ filename: string }>(
      'SELECT filename FROM schema_migrations',
    );
    return new Set(rows.map((r) => r.filename));
  }

  /** Returns the filenames that were applied by this run (empty if up to date). */
  async migrate(): Promise<string[]> {
    await this.ensureLedger();
    const applied = await this.appliedFilenames();

    const all = (await readdir(this.migrationsDir))
      .filter((f) => f.endsWith('.sql'))
      .sort();

    const pending = all.filter((f) => !applied.has(f));
    const executed: string[] = [];

    for (const filename of pending) {
      const sql = await readFile(path.join(this.migrationsDir, filename), 'utf8');
      await this.db.transaction(async (runner) => {
        await runner.query(sql);
        await runner.query('INSERT INTO schema_migrations (filename) VALUES ($1)', [filename]);
      });
      executed.push(filename);
    }

    return executed;
  }
}
