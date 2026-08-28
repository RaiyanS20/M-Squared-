# 01 — Architecture

> **Concept** → The code → Do it yourself → Check yourself

## The question every architecture answers

*When a requirement changes, how many files do I have to touch?*

That is all "good architecture" means. Not elegance, not patterns for their own
sake — just: does a small change stay small?

## Layers, and the one rule

The server has four layers. The rule that makes them worth having is about
**which way the arrows point**:

```
  ┌──────────────────────────────────────────────────────┐
  │  API layer        controllers, routing, error mapping │  knows HTTP
  └───────────────────────┬──────────────────────────────┘
                          │ calls
  ┌───────────────────────▼──────────────────────────────┐
  │  Service layer    TimelineService — the use cases     │  knows the product
  └───────────────────────┬──────────────────────────────┘
                          │ calls (through abstractions)
  ┌───────────────────────▼──────────────────────────────┐
  │  Repository layer CountryRepository (abstract)        │  knows storage
  │                   ├── PostgresCountryRepository       │
  │                   └── InMemoryCountryRepository       │
  └───────────────────────┬──────────────────────────────┘
                          │ uses
  ┌───────────────────────▼──────────────────────────────┐
  │  Domain layer     Year, Country, HistoricalEvent      │  knows the rules
  └──────────────────────────────────────────────────────┘
```

**Dependencies point inward. The domain depends on nothing.**

Verify it yourself — this is a rule you can check mechanically:

```bash
cd server
grep -rn "^import" src/domain/           # only imports its own siblings
grep -rn "express\|'pg'" src/services/   # no matches. The service knows neither.
grep -rln "from 'pg'" src/               # exactly one file: src/db/Database.js
```

`src/domain/Year.js` imports one thing: its sibling `errors.js`. It does not
know a database exists. That is not an accident; it is the whole design.

## Why bother? Three concrete payoffs

**1. The tests are fast because of it.** `TimelineService` accepts abstract
repositories, so tests hand it two arrays instead of a database. All 163 tests
run in about 1.5 seconds. A suite you can run on every save is a different tool
from one you run before lunch.

**2. Changing storage is a fifteen-line change.** Everything concrete is wired
in exactly one file, `src/main.js` — the *composition root*. Swapping PostgreSQL
for something else means writing a new repository class and editing those lines.
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
this course prepares you for, and it is why the structure is here despite the
app itself being small.

A useful way to hold it: *the number of layers should match how much the thing
is going to change.*

## Where the browser app fits

```
  index.html            the structure, written by hand once
        │
  TimelineView, CountryListView, EventListView   ← the only code touching the DOM
        │ told what to show by
  App.js                owns state, decides what to load
        │ calls
  TimelineApi           one method per feature
        │ extends
  ApiClient             fetch, parse, error mapping (one place)
        │ HTTP
  ══════════════════════ network boundary ══════════════════════
        │
  Express API
```

Same principle: **views never call `fetch`**, exactly as services never write
SQL. Chapters 7 and 8 cover this in detail.

## One server, one origin

Unusually for a modern project, there is no separate dev server:

```js
app.use('/api', new TimelineController(service).routes());
app.use(express.static(CLIENT_DIR));      // ← serves index.html, css/, js/
```

The browser loads the page and the API from the same origin, so there is no
proxy to configure and no CORS in development. Production splits them (chapter
11), which is where the CORS configuration finally matters.

## The request, end to end

Follow one click through the whole system. The user picks **1994**, then **South
Africa**:

| # | Where | What happens |
|---|-------|--------------|
| 1 | `CountryListView.js` | Click handler calls `onSelectCountry('ZAF')` |
| 2 | `App.js` | `store.setState({ selectedCountryCode: 'ZAF' })` |
| 3 | `Store.js` | Notifies subscribers; `App` starts the events request |
| 4 | `App.js` | Aborts any in-flight events request, sets status `loading` |
| 5 | `TimelineApi.getEvents` | `GET /api/years/1994/countries/ZAF/events` |
| 6 | `TimelineController` | Validates: is `1994` a real year? Is `ZAF` three letters? |
| 7 | `TimelineService` | Looks up the country, then its events. Missing → `NotFoundError` |
| 8 | `PostgresEventRepository` | Parameterised `SELECT ... WHERE year = $1 AND country_id = $2` |
| 9 | `RowMapper` | snake_case rows → `HistoricalEvent` domain objects |
| 10 | `TimelineController` | `res.json(...)` → each object's `toJSON()` shapes the wire format |
| 11 | `App.js` | State becomes `{ status: 'success', data }` |
| 12 | `EventListView` | Builds the cards with `document.createElement` |

**Twelve stops for one click. Is that not over-engineered?** Look at what each
stop is *for*: 4 makes a race condition impossible, 6 rejects bad input at the
edge, 7 is the only place product rules live, 8 makes SQL injection impossible,
9 stops database vocabulary leaking into the app, 10 stops private fields
leaking to the browser, 12 makes XSS impossible.

Each stop is one class of bug that cannot happen. In a single 200-line handler
doing all twelve jobs, every one of those is a thing you have to *remember*
every time.

## Do it yourself

1. **Trace it in reverse.** Pick the text "First non-racial democratic election"
   in the browser and find every file it passes through, from `seedData.js` to
   the screen. Write the list down.

2. **Break a layer on purpose.** Add `import pg from 'pg';` to
   `src/services/TimelineService.js`. Nothing stops you — JavaScript will not
   complain. Now ask what the project lost. (Run `npm test` and think about what
   would happen if that import were actually *used*.) Then remove it.

3. **Find the composition root.** Open `src/main.js` and list every concrete
   class named there. Then confirm no other file outside the wiring and the
   storage layer names one:
   ```bash
   cd server
   grep -rn "PostgresDatabase\|PostgresCountryRepository\|PostgresEventRepository" src/ \
     | grep -vE "src/main\.js|src/repositories/postgres/|src/db/|src/scripts/" \
     | grep -v '^\S*: *\*'
   ```
   You should get nothing. Note what the exclusions say: `src/scripts/` is
   excluded because `migrate.js` and `seed.js` are **their own composition
   roots** — separate entry points that wire their own small object graph. The
   last `grep -v` drops comment lines, which mention the class name in prose.

   That is worth pausing on: a project can have more than one composition root,
   one per entry point. What it must not have is concrete wiring scattered
   through the service and API layers.

## Check yourself

- Which layer would you change to support years before 100 AD?
- Which layer would you change to return XML instead of JSON?
- Which layer would you change to move from PostgreSQL to SQLite?
- Which layer would you change to require an API key?

<details>
<summary>Answers</summary>

- **Before 100 AD** — the domain only: `Year.EARLIEST` in `src/domain/Year.js`
  (plus the `events_year_in_range` CHECK constraint in a new migration, since
  the database enforces the same rule independently).
- **XML instead of JSON** — the API layer only. The controller decides the wire
  format; nothing below it knows or cares.
- **PostgreSQL → SQLite** — a new repository implementation, plus the wiring in
  `main.js`. The contract test suite runs against the new class unchanged and
  tells you whether you got it right.
- **API key** — the API layer, as Express middleware in `createApp`. Note this
  one is *cross-cutting*: it applies to every route rather than living inside any
  single one. That is what middleware is for.

If your instinct on all four was "the layer that owns that concern, and only
that layer", the architecture is doing its job.
</details>

→ Next: [02 — JavaScript and OOP](02-javascript-and-oop.md)
