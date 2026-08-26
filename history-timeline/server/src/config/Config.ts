import { z } from 'zod';

/**
 * Validates process.env once, at boot, and turns it into a typed object.
 *
 * The rule this enforces: a misconfigured server should fail LOUDLY on startup,
 * not mysteriously on the first request at 3am. Reading `process.env` scattered
 * through the codebase is how you end up with `undefined` in a connection string.
 */
const EnvSchema = z.object({
  NODE_ENV: z.enum(['development', 'test', 'production']).default('development'),
  PORT: z.coerce.number().int().positive().default(4000),
  DATABASE_URL: z.string().min(1, 'DATABASE_URL is required'),
  DATABASE_SSL: z
    .string()
    .optional()
    .transform((v) => v === 'true'),
  CORS_ORIGINS: z.string().default('http://localhost:5173'),
});

export class Config {
  private constructor(
    readonly nodeEnv: 'development' | 'test' | 'production',
    readonly port: number,
    readonly databaseUrl: string,
    readonly databaseSsl: boolean,
    readonly corsOrigins: string[],
  ) {}

  /**
   * OOP lesson: a PRIVATE CONSTRUCTOR plus a STATIC FACTORY. It is impossible to
   * construct a Config that skipped validation, because `new Config(...)` is not
   * reachable from outside this file.
   */
  static fromEnv(env: NodeJS.ProcessEnv = process.env): Config {
    const parsed = EnvSchema.safeParse(env);
    if (!parsed.success) {
      const issues = parsed.error.issues
        .map((i) => `  - ${i.path.join('.') || '(root)'}: ${i.message}`)
        .join('\n');
      throw new Error(`Invalid environment configuration:\n${issues}`);
    }
    const e = parsed.data;
    // In production the DB link crosses the network, so require TLS unless
    // explicitly told otherwise.
    const ssl = e.DATABASE_SSL || e.NODE_ENV === 'production';
    return new Config(
      e.NODE_ENV,
      e.PORT,
      e.DATABASE_URL,
      ssl,
      e.CORS_ORIGINS.split(',').map((o) => o.trim()).filter(Boolean),
    );
  }

  get isProduction(): boolean {
    return this.nodeEnv === 'production';
  }
}
