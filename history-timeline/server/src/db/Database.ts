import pg from 'pg';

/**
 * The narrow slice of "a database" that our repositories actually need.
 *
 * OOP lesson: PROGRAM TO AN INTERFACE. Repositories accept a `QueryRunner`, not
 * a `pg.Pool`. That keeps the `pg` dependency in exactly one file, and it means
 * a repository can be handed a pooled connection, a single transaction client,
 * or a fake — all three satisfy this interface.
 */
export interface QueryRunner {
  query<T extends pg.QueryResultRow = pg.QueryResultRow>(
    sql: string,
    params?: readonly unknown[],
  ): Promise<T[]>;
}

/**
 * "Something that can report whether it is usable."
 *
 * The health endpoint needs exactly one method, so this is what it depends on —
 * not the concrete `PostgresDatabase`. Interface segregation again: the API
 * layer stays free of any storage-specific type, and a test can pass `null` or
 * a two-line fake.
 */
export interface HealthCheck {
  isHealthy(): Promise<boolean>;
}

export interface DatabaseConfig {
  connectionString: string;
  /** Enable TLS. Required by AWS RDS in production. */
  ssl?: boolean;
  maxConnections?: number;
}

/**
 * Owns the PostgreSQL connection pool for the whole process.
 *
 * A pool is expensive to create and cheap to share, so exactly one of these is
 * constructed in `main.ts` and passed down. That is dependency injection done by
 * hand — no framework, no magic, and completely explicit about who owns what.
 */
export class PostgresDatabase implements QueryRunner, HealthCheck {
  private readonly pool: pg.Pool;

  constructor(config: DatabaseConfig) {
    this.pool = new pg.Pool({
      connectionString: config.connectionString,
      max: config.maxConnections ?? 10,
      // RDS presents a certificate signed by the Amazon RDS CA. For a learning
      // project we accept it without pinning the CA bundle; the deployment
      // chapter (docs/10) explains how to do this properly in production.
      ssl: config.ssl ? { rejectUnauthorized: false } : false,
      // Fail fast instead of hanging a request forever if the DB is unreachable.
      connectionTimeoutMillis: 5_000,
    });

    // An idle client erroring out (e.g. RDS failover) must not crash the process.
    this.pool.on('error', (err: Error) => {
      console.error('[db] idle client error', err);
    });
  }

  async query<T extends pg.QueryResultRow = pg.QueryResultRow>(
    sql: string,
    params: readonly unknown[] = [],
  ): Promise<T[]> {
    const result = await this.pool.query<T>(sql, params as unknown[]);
    return result.rows;
  }

  /**
   * Runs `work` inside a transaction on a single dedicated client, rolling back
   * if anything throws. Used by the migrator and the seed script.
   */
  async transaction<T>(work: (runner: QueryRunner) => Promise<T>): Promise<T> {
    const client = await this.pool.connect();
    const runner: QueryRunner = {
      query: async (sql, params = []) => (await client.query(sql, params as unknown[])).rows as never,
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
      // Always return the client to the pool, success or failure.
      client.release();
    }
  }

  /** Used by the /health endpoint and by the AWS load balancer health check. */
  async isHealthy(): Promise<boolean> {
    try {
      await this.query('SELECT 1');
      return true;
    } catch {
      return false;
    }
  }

  async close(): Promise<void> {
    await this.pool.end();
  }
}
