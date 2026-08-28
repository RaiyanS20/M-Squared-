# 05 — The repository pattern

> **Concept** → The code → Do it yourself → Check yourself

## The problem it solves

Without a repository, SQL leaks upward:

```js
// The route now knows about SQL, snake_case, and connection pools.
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
   browser — the database schema *is* the public API, so every column rename is
   a breaking change.
3. **The next endpoint copies the query** and the two drift apart.

## The pattern

A **repository** is an object that looks like a collection of domain objects and
hides where they are actually stored.

```js
export class CountryRepository {
  async findAll() { ... }
  async findAllWithEventsInYear(year) { ... }
  async findById(id) { ... }
  async findByCode(code) { ... }
}
```

No SQL. No `pg`. No HTTP. It returns **domain objects**, not rows — that is the
part people skip, and it is the part that matters. A repository returning raw
rows has moved the problem, not solved it.

## Two implementations, one contract

```
              CountryRepository (abstract)
                      ▲
        ┌─────────────┴─────────────┐
PostgresCountryRepository    InMemoryCountryRepository
   real SQL, real pool          arrays, used by tests
```

`TimelineService` accepts either. That single fact is why the service tests run
in milliseconds.

## `RowMapper`: the anti-corruption layer

```js
static toEvent(row) {
  return new HistoricalEvent({
    id: row.id,
    countryId: row.country_id,          // snake_case → camelCase
    year: new Year(row.year),           // number → validated value object
    monthDay: row.month_day?.trim() || null,   // CHAR(5) is blank-padded
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

Here is what makes this pattern trustworthy rather than merely tidy.

Two implementations of one abstraction is a **promise**: anywhere a
`CountryRepository` is expected, either works. Promises rot. One sorts by name,
the other returns insertion order; tests pass, production is subtly wrong.

**This matters more in JavaScript than in a typed language.** A compiler would at
least tell you a method was missing. Here, nothing checks that the two agree
until something calls them — so the tests have to be that something.

`tests/repositoryContract.test.js`:

```js
function contractFor(label, build, expectations) {
  describe(`${label} satisfies the repository contract`, () => {
    it('findAll returns countries sorted by name', ...);
    it('findByCode is case-insensitive', ...);
    it('findByCode returns null for an unknown code', ...);
    // ...
  });
}

contractFor('InMemory repositories', ...);                    // always runs
describe('PostgreSQL repositories', { skip: !testDatabaseUrl }, () => {
  contractFor('Postgres repositories', ...);                  // opt-in
});
```

Nine assertions × two implementations. The in-memory pair always runs, so the
fast suite stays fast. The Postgres pair runs when you ask:

```bash
npm run db:up && npm run db:migrate && npm run db:seed
TEST_DATABASE_URL=postgres://historian:historian@localhost:5432/history_timeline npm test
```

This generalises: **whenever you have an interface with more than one
implementation, write the tests against the interface.** It is the difference
between "we have an abstraction" and "our abstraction is true".

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

**Returning `null`, not throwing.** `findByCode` returns a `Country` or `null`.
"Not found" is a normal outcome for a lookup; it is the *service* that decides
whether that is an error (`NotFoundError`) or fine. Repositories report;
services judge.

**Copy before sorting.** `InMemoryCountryRepository.findAll` does
`[...this.#data.countries].sort(...)`, because `Array.prototype.sort` **mutates
in place** and that array belongs to the dataset. A "read" method that quietly
reorders shared data is a genuinely confusing bug.

## One more inversion

`PostgresCountryRepository` takes a `db` — anything with a `query(sql, params)`
method — not a `pg.Pool`. So it can be handed a pooled connection, a
transaction's client, or a fake, and the `pg` import lives in exactly one file:

```bash
grep -rln "from 'pg'" server/src/     # one file
```

The transaction helper relies on this: `db.transaction()` hands the callback a
`runner` with the same `query` shape, so the seed script's writes go through the
transaction's single client without any code knowing the difference.

## Do it yourself

1. **Break the contract, watch it caught.** Delete `.sort(...)` from
   `InMemoryCountryRepository.findAll()`. Run
   `node --test server/tests/repositoryContract.test.js`. Restore it.

2. **Add a method to the contract.** Add `findByRegion(region)` to
   `CountryRepository` (throwing the "must implement" error), then implement it
   in both. Add contract assertions and confirm both pass. Notice the base class
   telling you exactly what is missing when you forget one.

3. **Write a third implementation.** Create a `CachingCountryRepository` that
   wraps another repository and memoises `findAll()`. Run it through
   `contractFor` — if the contract passes, it is safe to substitute anywhere.
   (This is the *decorator* pattern, and it composes precisely because everything
   agrees on one contract.)

4. **Cause an N+1 on purpose.** Rewrite `getTimelineScale` to call
   `findByYear(y)` for each year in the range. Time both with `console.time`.
   Then put it back.

## Check yourself

- Why do repositories return domain objects rather than rows?
- Why is `findByCode` case-insensitive at the repository level, not in the
  controller?
- When is the repository pattern *not* worth it?

<details>
<summary>Answers</summary>

- **Domain objects**: a row is unvalidated data in the database's vocabulary.
  Returning rows pushes mapping and validation onto every caller and re-couples
  the layers you just separated.
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
