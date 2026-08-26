# 01 — Architecture

> **Concept** → The code → Do it yourself → Check yourself

## The question every architecture answers

*When a requirement changes, how many files do I have to touch?*

That is all "good architecture" means. Not elegance, not patterns for their own
sake — just: does a small change stay small?

## Layers, and the one rule

This project has four layers on the server. The rule that makes them worth
having is about **which way the arrows point**:

```
  ┌──────────────────────────────────────────────────────┐
  │  API layer        controllers, routing, error mapping│  knows HTTP
  └───────────────────────┬──────────────────────────────┘
                          │ calls
  ┌───────────────────────▼──────────────────────────────┐
  │  Service layer    TimelineService — the use cases     │  knows the product
  └───────────────────────┬──────────────────────────────┘
                          │ calls (through abstractions)
  ┌───────────────────────▼──────────────────────────────┐
  │  Repository layer abstract CountryRepository etc.     │  knows storage
  │                   ├── PostgresCountryRepository       │
  │                   └── InMemoryCountryRepository       │
  └───────────────────────┬──────────────────────────────┘
                          │ uses
  ┌───────────────────────▼──────────────────────────────┐
  │  Domain layer     Year, Country, HistoricalEvent      │  knows the rules
  └──────────────────────────────────────────────────────┘
```

**Dependencies point inward. The domain depends on nothing.**

Check it yourself — this is a rule you can verify mechanically:

```bash
cd server
grep -rn "import" src/domain/          # only imports other domain files
grep -rn "express\|pg" src/services/   # no matches. The service knows neither.
```

`src/domain/Year.ts` has no imports at all except its sibling `errors.ts`. It
does not know a database exists. That is not an accident; it is the whole design.

## Why bother? Three concrete payoffs

**1. The tests are fast because of it.** `TimelineService` accepts abstract
repositories, so tests hand it two arrays instead of a database. All 97 server
tests run in about one second. A test suite you can run on every save is a
different tool from one you run before lunch.

**2. Changing storage is a fifteen-line change.** Everything concrete is wired
in exactly one file, `src/main.ts` — the *composition root*. Swapping PostgreSQL
for DynamoDB means writing a new repository class and editing those lines.
Nothing in `services/`, `domain/` or `api/` changes.

**3. A rule lives in exactly one place.** "An event must cite a source" is in
`HistoricalEvent`'s constructor. The seed script, the API and every test get it
for free, because none of them can build an event any other way.

## The trade-off, stated honestly

This structure has a real cost: **more files, and more indirection**. Reading
"list the countries for a year" means opening a controller, a service, an
abstract repository and a concrete repository. For a weekend project that ships
once and is never touched again, that is genuine overhead and you should not do
it.

It pays off when a project is **maintained** — when requirements change, when
more than one person works on it, when it has to be tested. That is the case
this course is preparing you for, and it is why the structure is here despite
the app itself being small.

A useful way to hold it: *the number of layers should match how much the thing
is going to change.* This one is deliberately at the "industry-standard,
maintained service" end so you see the whole pattern.

## Where the frontend fits

The client is a separate program that talks to the server over JSON:

```
  React components      ← what the user sees
        │ call
  hooks (useTimeline)   ← loading / error / data states
        │ call
  TimelineApi           ← one method per feature
        │ extends
  ApiClient             ← fetch, parse, error mapping (one place)
        │ HTTP
  ══════════════════════════ network boundary ══════════════════════
        │
  Express API
```

Same principle: components never call `fetch`, exactly as services never write
SQL. Chapter 7 covers this in detail.

## The request, end to end

Follow one click through the whole system. The user picks **1994**, then **South
Africa**:

| # | Where | What happens |
|---|-------|--------------|
| 1 | `CountryPicker.tsx` | Button click calls `onSelectCountry('ZAF')` |
| 2 | `App.tsx` | `setSelectedCountryCode('ZAF')` — state changes, React re-renders |
| 3 | `useEvents` | Dependencies changed → aborts any in-flight request, starts a new one |
| 4 | `TimelineApi.getEvents` | `GET /api/years/1994/countries/ZAF/events` |
| 5 | `TimelineController` | Validates: is `1994` a real year? Is `ZAF` three letters? |
| 6 | `TimelineService` | Looks up the country, then its events. Country missing → `NotFoundError` |
| 7 | `PostgresEventRepository` | Parameterised `SELECT ... WHERE year = $1 AND country_id = $2` |
| 8 | `RowMapper` | snake_case rows → `HistoricalEvent` domain objects |
| 9 | `TimelineController` | `res.json(...)` → each object's `toJSON()` shapes the wire format |
| 10 | `useEvents` | State becomes `{ status: 'success', data }` |
| 11 | `EventList.tsx` | Renders the cards |

**Ten stops for one click. Is that not over-engineered?** Look at what each stop
is *for*: 5 rejects bad input at the edge, 6 is the only place product rules
live, 7 makes SQL injection impossible, 8 stops database vocabulary leaking into
the app, 9 stops private fields leaking to the browser. Each stop is one class of
bug that cannot happen. In a single 200-line handler doing all ten jobs, every
one of those is a thing you have to *remember* every time.

## Do it yourself

1. **Trace it in reverse.** Pick the text "First non-racial democratic election"
   in the browser and find every file it passes through, from `seedData.ts` to
   the screen. Write the list down.

2. **Break a layer on purpose.** Add `import pg from 'pg';` to
   `src/services/TimelineService.ts`. Nothing stops you — the compiler is fine
   with it. Now ask: what did the project lose? (Run `npx vitest run
   tests/unit/TimelineService.test.ts` and think about what would happen if that
   import were actually *used*.) Then remove it.

3. **Find the composition root.** Open `src/main.ts` and list every concrete
   class named there. Then confirm that no file outside the wiring and the
   storage layer names a concrete `Postgres*` type:
   ```bash
   cd server
   grep -rn "PostgresDatabase\|PostgresCountryRepository\|PostgresEventRepository" src/ --include=*.ts \
     | grep -vE "src/main\.ts|src/repositories/postgres/|src/db/"
   ```
   You should get nothing. The API layer depends on a two-line `HealthCheck`
   interface rather than on `PostgresDatabase` for exactly this reason — an
   earlier draft of this project failed that check, and the grep is how it was
   found.

## Check yourself

- Which layer would you change to support years before 100 AD?
- Which layer would you change to return XML instead of JSON?
- Which layer would you change to move from PostgreSQL to SQLite?
- Which layer would you change to require an API key?

<details>
<summary>Answers</summary>

- **Before 100 AD** — the domain only: `Year.EARLIEST` in `src/domain/Year.ts`
  (plus the `events_year_in_range` CHECK constraint in a new migration, since the
  database enforces the same rule independently).
- **XML instead of JSON** — the API layer only. The controller decides the wire
  format; nothing below it knows or cares.
- **PostgreSQL → SQLite** — a new repository implementation, plus the wiring in
  `main.ts`. The contract test suite in `tests/unit/repositoryContract.test.ts`
  runs against the new class unchanged and tells you whether you got it right.
- **API key** — the API layer, as Express middleware in `createApp`. Note this
  one is *cross-cutting*: it applies to every route rather than living inside
  any single one. That is what middleware is for.

If your instinct on all four was "the layer that owns that concern, and only
that layer", the architecture is doing its job.
</details>

→ Next: [02 — TypeScript and OOP](02-typescript-and-oop.md)
