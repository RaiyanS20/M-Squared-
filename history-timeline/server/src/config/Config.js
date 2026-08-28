/**
 * Validates the environment once, at boot, and returns a frozen object.
 *
 * The rule this enforces: a misconfigured server should fail LOUDLY on startup,
 * not mysteriously on the first request at 3am. Reading `process.env` scattered
 * through a codebase is how you end up with the string "undefined" inside a
 * database connection string.
 *
 * OOP lesson: PRIVATE CONSTRUCTOR + STATIC FACTORY. `new Config()` throws, so
 * `Config.fromEnv()` is the only way in — which means a Config that skipped
 * validation cannot exist anywhere in the program.
 *
 * (A typed project might reach for a schema library here. In plain JavaScript,
 * writing the twenty lines yourself is clearer and adds no dependency.)
 */
export class Config {
  static #buildToken = Symbol('Config.build');

  constructor(token, values) {
    if (token !== Config.#buildToken) {
      throw new TypeError('Use Config.fromEnv() — the constructor is private.');
    }
    Object.assign(this, values);
    Object.freeze(this);
  }

  static fromEnv(env = process.env) {
    const problems = [];

    const nodeEnv = env.NODE_ENV ?? 'development';
    if (!['development', 'test', 'production'].includes(nodeEnv)) {
      problems.push(`NODE_ENV must be development, test or production (got "${nodeEnv}")`);
    }

    const port = Number.parseInt(env.PORT ?? '4000', 10);
    if (!Number.isInteger(port) || port < 1 || port > 65535) {
      problems.push(`PORT must be a number between 1 and 65535 (got "${env.PORT}")`);
    }

    const databaseUrl = env.DATABASE_URL ?? '';
    if (!databaseUrl) {
      problems.push('DATABASE_URL is required (see .env.example)');
    } else if (!/^postgres(ql)?:\/\//.test(databaseUrl)) {
      problems.push('DATABASE_URL must start with postgres:// or postgresql://');
    }

    if (problems.length > 0) {
      throw new Error(
        `Invalid environment configuration:\n${problems.map((p) => `  - ${p}`).join('\n')}`,
      );
    }

    return new Config(Config.#buildToken, {
      nodeEnv,
      port,
      databaseUrl,
      // In production the DB link crosses a network, so require TLS unless
      // explicitly told otherwise.
      databaseSsl: env.DATABASE_SSL === 'true' || nodeEnv === 'production',
      corsOrigins: (env.CORS_ORIGINS ?? '')
        .split(',')
        .map((o) => o.trim())
        .filter(Boolean),
      isProduction: nodeEnv === 'production',
    });
  }
}
