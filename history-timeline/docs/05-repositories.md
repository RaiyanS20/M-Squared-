# 05 — The repository pattern

> **Concept** → The code → Do it yourself → Check yourself

## The problem it solves

Without a repository, SQL leaks upward:

```ts
// The controller now knows about SQL, snake_case, and connection pools.
app.get('/api/years/:year/countries/:code/events', async (req, res) => {
  const rows = await pool.query(
    `SELECT e.* FROM historical_events e
       JOIN countries c ON c.id = e.country_id
      WHERE e.year = $1 AND c.code = $2`, [req.params.year, req.params.code]);
  res.json(rows.rows);
});
```

Three problems, and they compound:

1. **Untestable** without a real database.
2. **`res.json(rows.rows)`** ships `country_id` and `created_at` straight to the
   browser — the database schema *is* the public API, so every column rename is a
   breaking change.
3. **The next endpoint copies the query** and the two drift apart.

## The pattern

A **repository** is an object that looks like a collection of domain objects and
hides where they are actually stored.

```ts
export abstract class CountryRepository {
  abstract findAll(): Promise<Country[]>;
  abstract findAllWithEventsInYear(year: number): Promise<Country[]>;
  abstract findById(id: number): Promise<Country | null>;
  abstract findByCode(code: string): Promise<Country | null>;
}
```

No SQL. No `pg`. No HTTP. It returns **domain objects**, not rows — that is the
part people skip, and it is the part that matters. A repository that returns raw
rows has moved the problem, not solved it.

## Two implementations, one contract

```
              CountryRepository (abstract)
                      ▲
        ┌─────────────┴─────────────┐
PostgresCountryRepository    InMemoryCountryRepository
   real SQL, real pool          arrays, used by tests
```

`TimelineService` accepts the abstract type, so it works with either. That single
fact is why the service tests run in milliseconds.

## `RowMapper`: the anti-corruption layer

```ts
static toEvent(row: EventRow): HistoricalEvent {
  return new HistoricalEvent({
    id: row.id,
    countryId: row.country_id,        // snake_case → camelCase
    year: new Year(row.year),         // number → validated value object
    category: row.category as EventCategory,
    monthDay: row.month_day?.trim() || null,  // CHAR(5) is blank-padded
    ...
  });
}
```

The database speaks `snake_case` and knows nothing about `Year`. The domain
speaks `camelCase` and refuses invalid data. **This is the only place the two
vocabularies meet**, so renaming a column touches one file.

Note `row.code.trim()` in `toCountry`. `CHAR(3)` is blank-padded by PostgreSQL,
so `'FRA'` can come back as `'FRA '`. Find that bug once at the boundary, fix it
once, and it never bites again — that is what a mapping layer is for.

## Contract tests: making Liskov real

Here is the thing that makes this pattern trustworthy rather than merely tidy.

Two implementations of one abstraction is a **promise**: anywhere a
`CountryRepository` is expected, either one works. Promises rot. One sorts by
name, the other returns insertion order; tests pass, production is subtly wrong.

So the assertions live in one function, and both implementations are fed through
it — `tests/unit/repositoryContract.test.ts`:

```ts
function contractFor(label, build, expectations) {
  describe(`${label} satisfies the repository contract`, () => {
    it('findAll returns countries sorted by name', async () => { ... });
    it('findByCode is case-insensitive', async () => { ... });
    it('findByCode returns null for an unknown code', async () => { ... });
    // ...
  });
}

contractFor('InMemory repositories', ...);   // always runs
describe.skipIf(!testDatabaseUrl)('PostgreSQL repositories', () => {
  contractFor('Postgres repositories', ...); // opt-in
});
```

Nine assertions × two implementations. The in-memory pair always runs, so the
fast suite stays fast. The Postgres pair runs when you ask for it:

```bash
npm run db:up && npm run db:migrate && npm run db:seed
TEST_DATABASE_URL=postgres://historian:historian@localhost:5432/history_timeline npm test
```

This is the pattern's real payoff, and it generalises: **whenever you have an
interface with more than one implementation, write the tests against the
interface.** It is the difference between "we have an abstraction" and "our
abstraction is true".

## Design notes on the specific methods

**`findAllWithEventsInYear`** exists because a list of 195 countries where 190
lead to an empty page is a bad experience. It uses `EXISTS` rather than a `JOIN`
+ `DISTINCT`:

```sql
SELECT c.id, c.code, c.name, c.region FROM countries c
 WHERE EXISTS (SELECT 1 FROM historical_events e
                WHERE e.country_id = c.id AND e.year = $1)
```

`EXISTS` stops at the first match per country and cannot produce duplicates, so
there is no `DISTINCT` to deduplicate afterwards.

**`countByYearRange`** returns counts for a whole range in **one query**. The
naive alternative — one query per year — is 126 round trips to render one page.

```
One query returning 60 rows:   ~2 ms
126 queries returning 1 row:   ~250 ms
```

This is the **N+1 query problem**, the most common performance bug in
database-backed applications, and it is invisible on a laptop with a local
database. Watch for it any time you see a query inside a loop.

**Returning `null`, not throwing.** `findByCode` returns `Country | null`. "Not
found" is a normal outcome for a lookup; it is the *service* that decides
whether that is an error (`NotFoundError`) or fine. Repositories report; services
judge.

## `QueryRunner`: one more inversion

```ts
export interface QueryRunner {
  query<T>(sql: string, params?: readonly unknown[]): Promise<T[]>;
}
```

The Postgres repositories accept a `QueryRunner`, not a `pg.Pool`. So they can be
handed a pooled connection, a single transaction client, or a fake — and the
`pg` import lives in exactly one file (`db/Database.ts`). Verify:

```bash
grep -rln "from 'pg'" server/src/    # one file
```

## Do it yourself

1. **Break the contract, watch it caught.** Delete `.sort(...)` from
   `InMemoryCountryRepository.findAll()`. Run `npx vitest run repositoryContract`.
   Restore it.

2. **Add a method to the contract.** Add `findByRegion(region: string)` to
   `CountryRepository`. The compiler will immediately fail both implementations
   until you write them — that is the abstract class earning its keep. Add
   contract assertions for it and confirm both pass.

3. **Write a third implementation.** Create a `CachingCountryRepository` that
   wraps another repository and memoises `findAll()`. Run it through
   `contractFor` — if the contract passes, it is safe to substitute anywhere.
   (This is the *decorator* pattern, and it composes precisely because everything
   agrees on one contract.)

4. **Cause an N+1 on purpose.** Rewrite `getTimelineScale` to call
   `findByYear(y)` for each year in the range. Time both versions with
   `console.time`. Then put it back.

## Check yourself

- Why do repositories return domain objects rather than rows?
- Why is `findByCode` case-insensitive at the repository level and not in the
  controller?
- When is the repository pattern *not* worth it?

<details>
<summary>Answers</summary>

- **Domain objects**: because a row is unvalidated data in the database's
  vocabulary. Returning rows pushes mapping and validation onto every caller,
  and re-couples the layers you just separated.
- **Case-insensitivity belongs in the repository**: it is a property of *how
  country codes are looked up*, true for every caller, not just for HTTP. Put it
  in the controller and the next caller (a CLI, a background job) gets different
  behaviour. Both implementations are tested for it in the contract.
- **Not worth it**: a script, a prototype, or anything where the data access will
  never be tested, swapped or reused. The pattern buys testability and
  substitutability; if you need neither, it is ceremony. Be honest about which
  situation you are in.

</details>

→ Next: [06 — The HTTP API](06-http-api.md)
