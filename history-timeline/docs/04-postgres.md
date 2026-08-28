# 04 — PostgreSQL

> **Concept** → The code → Do it yourself → Check yourself

## Why a relational database here

Our data is **relational**: events belong to countries, and the core question is
"events where year = X and country = Y". PostgreSQL gives us three things that
matter for a factual history project:

- **Constraints** — the database itself refuses bad data.
- **Transactions** — the seed either fully applies or does not apply at all.
- **Real query planning** — the index makes the hot query fast, and `EXPLAIN`
  proves it rather than you guessing.

## The schema

`server/src/db/migrations/001_create_countries_and_events.sql`. Every line is a
decision:

```sql
CREATE TABLE countries (
    id      SERIAL PRIMARY KEY,
    code    CHAR(3)      NOT NULL UNIQUE,
    name    TEXT         NOT NULL UNIQUE,
    region  TEXT         NOT NULL,
    CONSTRAINT countries_code_is_alpha CHECK (code ~ '^[A-Z]{3}$')
);
```

**Why a surrogate `id` when `code` is already unique?** `code` is a *natural key*
— it has real-world meaning, which is exactly why it can change (countries are
renamed, ISO reassigns codes). A meaningless integer is stable forever, keeps
foreign keys narrow, and means renaming a country does not rewrite every event.

**Why `year INTEGER` rather than `DATE`?** As chapter 3 covered: `DATE` cannot
represent "sometime in 1347" and would force us to invent a day.

```sql
country_id INTEGER NOT NULL REFERENCES countries (id) ON DELETE CASCADE
```

**`ON DELETE CASCADE`**: delete a country, its events go too. The alternative is
orphan rows that no query reaches but every `COUNT(*)` includes. Choose the
behaviour deliberately — the other options are `RESTRICT` (refuse the delete) and
`SET NULL`.

## Constraints: the last line of defence

```sql
CONSTRAINT events_source_is_url  CHECK (source_url ~ '^https?://'),
CONSTRAINT events_category_known CHECK (category IN ('politics', 'conflict', ...)),
CONSTRAINT events_unique_per_country_year UNIQUE (country_id, year, title)
```

These duplicate rules already in `HistoricalEvent`. **That duplication is
correct.** Application code can be bypassed; the database cannot be bypassed by
anything that writes to the database.

Try it:

```bash
psql history_timeline -c \
  "INSERT INTO historical_events (country_id, year, title, summary, category, source_url)
   VALUES (1, 1969, 'Made up', 'No source', 'politics', 'not-a-url');"
-- ERROR:  new row for relation "historical_events" violates check constraint "events_source_is_url"
```

The editorial promise holds even against someone with a database password.

## Indexes: make the hot query fast

```sql
CREATE INDEX idx_events_year_country ON historical_events (year, country_id);
```

Without an index, finding events for 1994 means reading **every row**. With 124
rows that is instant; with 10 million it is a timeout.

**Column order matters.** An index on `(year, country_id)` serves:

- ✅ `WHERE year = 1994 AND country_id = 12` — both columns
- ✅ `WHERE year = 1994` — leading column alone
- ❌ `WHERE country_id = 12` — *not* the leading column

Think of a phone book sorted by (surname, first name): useless for finding
everyone called "James".

Prove it on your own machine:

```sql
EXPLAIN ANALYZE SELECT * FROM historical_events WHERE year = 1994 AND country_id = 12;
```

With 124 rows PostgreSQL may still choose a sequential scan — reading a tiny
table is cheaper than consulting an index. That is the planner being smart, not
the index being wrong. To see the index chosen, generate volume:

```sql
INSERT INTO historical_events (country_id, year, title, summary, category, source_url)
SELECT 1, 1900 + (n % 120), 'Filler ' || n, 'Generated row for an index demo.',
       'politics', 'https://example.org/' || n
FROM generate_series(1, 200000) AS n;

ANALYZE historical_events;
EXPLAIN ANALYZE SELECT * FROM historical_events WHERE year = 1994 AND country_id = 1;
-- now: Index Scan using idx_events_year_country
```

Clean up with `DELETE FROM historical_events WHERE title LIKE 'Filler %';`.

**Measure, do not guess.** Indexes cost write speed and disk. Add them for
queries you actually run.

## Migrations

`src/db/Migrator.js` — about thirty lines. A ledger table, a sorted file list, a
transaction per file:

```sql
CREATE TABLE IF NOT EXISTS schema_migrations (filename TEXT PRIMARY KEY, ...)
```

The rules that keep a team out of trouble:

1. **Never edit an applied migration.** Once `001_*.sql` has run anywhere but
   your laptop, it is history. Fix it with `002_*.sql`.
2. **One transaction per file**, so a failure leaves no half-applied state.
3. **Zero-pad filenames** (`001_`, `002_`) so lexical order is real order.
4. **Migrations run before the new code deploys.** Chapter 11 covers sequencing.

Writing this by hand is deliberate: every migration tool you meet later (Flyway,
Prisma Migrate, Knex) is this plus features, and you now know what the features
are *for*.

## SQL injection, and the one habit that prevents it

```js
// NEVER. This is how databases get dumped.
db.query(`SELECT * FROM countries WHERE code = '${code}'`);
//  code = "' OR '1'='1"                       → returns everything
//  code = "'; DROP TABLE historical_events; --" → exactly what it looks like

// ALWAYS. The value can never be parsed as SQL.
db.query('SELECT * FROM countries WHERE code = $1', [code]);
```

`pg` sends the query and the parameters to PostgreSQL **separately**. The server
plans the SQL first, then binds values — a value is data, and there is no
mechanism by which it becomes code.

Audit it yourself — look for any interpolation in a SQL string:

```bash
cd server && grep -n '\${' src/repositories/postgres/*.js
```

Two hits, both `${EVENT_COLUMNS}` — a hard-coded constant at the top of the same
file, never user input. Every actual *value* goes through `$1`, `$2`.

That grep is worth keeping: **"is there interpolation in any SQL string?" is a
review question with a yes/no answer**, and it catches injection before it ships.

## Connection pooling

`PostgresDatabase` wraps a `pg.Pool`. Opening a PostgreSQL connection costs
several milliseconds and real server memory; a pool keeps a handful open and
lends them out.

```js
max: 10,                        // connections in this process
connectionTimeoutMillis: 5000,  // fail fast rather than hang forever
```

**Sizing matters on AWS.** `max` is per *process*. Four containers × 10 = 40
connections. An RDS `db.t4g.micro` allows roughly 80. Scale to twenty without
thinking and you exhaust connections before CPU — a genuinely confusing outage.
Chapter 11 returns to this.

The `pool.on('error')` handler matters too: an idle client erroring (an RDS
failover) would otherwise crash the whole process.

## Transactions

```js
await db.transaction(async (runner) => {
  for (const country of SEED_COUNTRIES) { await runner.query(...); }
  for (const event of SEED_EVENTS)     { await runner.query(...); }
});
```

All of it, or none of it. Note the implementation detail in `Database.js`: the
transaction runs on **one dedicated client** from the pool (`pool.connect()`),
not on the pool itself. Issuing `BEGIN` on a pool would start a transaction on
one random connection and the next statement might land on a different one.
`finally { client.release() }` returns it whatever happens — forget that and the
pool leaks connections until the app stops responding.

## Idempotent seeds

```sql
INSERT INTO countries (code, name, region) VALUES ($1, $2, $3)
ON CONFLICT (code) DO UPDATE SET name = EXCLUDED.name, region = EXCLUDED.region
```

Run the seed twice, get the same 124 rows — not 248. So you can correct a typo
in `seedData.js`, re-run, and the fix lands. **Any script you might run twice
should be safe to run twice**; it is one of the cheapest reliability habits there
is.

## Do it yourself

1. **Watch a constraint fire.** Run the bad `INSERT` above. Then try a duplicate
   `(country_id, year, title)`. Then a `category` of `'gossip'`.

2. **Prove the seed is idempotent.**
   ```bash
   npm run db:seed && npm run db:seed
   psql history_timeline -c 'SELECT COUNT(*) FROM historical_events;'   # 124
   ```

3. **Write a migration.** Add an optional `place TEXT` column in
   `002_add_event_place.sql`. Run `npm run db:migrate`. Run it again — nothing
   should happen. Then check the ledger:
   `SELECT * FROM schema_migrations;`

4. **Ask the database a question it is good at.** Which decade has the most
   recorded events?
   ```sql
   SELECT (year / 10) * 10 AS decade, COUNT(*) AS events
     FROM historical_events GROUP BY decade ORDER BY events DESC LIMIT 5;
   ```

## Check yourself

- Why does `COUNT(*)` come back as a string from `pg`?
- Why is `month_day` `CHAR(5)` and not two integer columns?
- What happens to events when their country is deleted, and why is that right
  here?

<details>
<summary>Answers</summary>

- **`COUNT(*)` is a string**: PostgreSQL returns `BIGINT` (64-bit), which exceeds
  JavaScript's safe integer range. `pg` hands you a string rather than silently
  losing precision, which is why `PostgresEventRepository` calls
  `Number.parseInt` explicitly. A library choosing correctness over convenience.
- **`CHAR(5)`**: it is a single indivisible fact, it sorts correctly as text
  (`"03-02" < "04-28"`), and it maps straight to the domain's `monthDay` with no
  assembly. Two integers would need a `CHECK` on each plus reassembly everywhere.
  (The cost: `CHAR` is blank-padded, which is why `RowMapper` calls `.trim()`.)
- **`ON DELETE CASCADE`**: the events disappear. Right here because an event
  without a country is meaningless in this model. If events could meaningfully
  outlive a country you would want `RESTRICT` instead.

</details>

→ Next: [05 — The repository pattern](05-repositories.md)
