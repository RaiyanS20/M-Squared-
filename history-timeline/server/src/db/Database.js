import pg from 'pg';

/**
 * Owns the PostgreSQL connection pool for the whole process.
 *
 * A pool is expensive to create and cheap to share, so exactly ONE of these is
 * built in `main.js` and passed down. That is dependency injection done by hand
 * — no framework, no magic, completely explicit about who owns what.
 *
 * This is also the only file in the project that imports `pg`. Check it:
 *     grep -rln "from 'pg'" server/src/
 */
export class PostgresDatabase {
  #pool;

  constructor({ connectionString, ssl = false, maxConnections = 10 }) {
    this.#pool = new pg.Pool({
      connectionString,
      max: maxConnections,
      // RDS presents a certificate signed by the Amazon RDS CA. For a learning
      // project we accept it without pinning the CA bundle; docs/10 explains how
      // to do this properly in production.
      ssl: ssl ? { rejectUnauthorized: false } : false,
      // Fail fast instead of hanging a request forever if the DB is unreachable.
      connectionTimeoutMillis: 5000,
    });

    // An idle client erroring out (an RDS failover, say) must not crash the
    // whole process. Without this handler, it would.
    this.#pool.on('error', (err) => {
      console.error('[db] idle client error', err);
    });
  }

  /**
   * Runs a query and returns just the rows.
   *
   * `params` is separate from `sql` on purpose — see docs/04. This is what makes
   * SQL injection impossible rather than merely unlikely.
   */
  async query(sql, params = []) {
    const result = await this.#pool.query(sql, params);
    return result.rows;
  }

  /**
   * Runs `work` inside a transaction on a single dedicated client, rolling back
   * if anything throws.
   *
   * Note it takes a client from the pool rather than using the pool directly.
   * Issuing BEGIN on a pool would start a transaction on one random connection,
   * and the next statement might land on a different one.
   */
  async transaction(work) {
    const client = await this.#pool.connect();
    const runner = {
      query: async (sql, params = []) => (await client.query(sql, params)).rows,
    };
    try {
      await client.query('BEGIN');
      const result = await work(runner);
      await client.query('COMMIT');
      return result;
    } catch (error) {
      await client.query('ROLLBACK');
      throw error;
    } finally {
      // Always return the client to the pool, success or failure. Forget this
      // and the pool leaks connections until the app stops responding.
      client.release();
    }
  }

  /** Used by /ready and by the AWS load balancer health check. */
  async isHealthy() {
    try {
      await this.query('SELECT 1');
      return true;
    } catch {
      return false;
    }
  }

  async close() {
    await this.#pool.end();
  }
}
