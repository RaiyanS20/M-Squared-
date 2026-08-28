# History Timeline

A factual, educational history site — **and a project-based course in building
one with vanilla JavaScript, HTML and CSS.**

Pick a year on the timeline, choose a country, and read what is documented to
have happened there. Every event links to its source.

![JavaScript](https://img.shields.io/badge/JavaScript-ES2022-f7df1e) ![Node](https://img.shields.io/badge/Node-22-339933) ![Postgres](https://img.shields.io/badge/PostgreSQL-16-336791) ![No build step](https://img.shields.io/badge/build_step-none-brightgreen)

---

## Two things at once

**A working application.** Vanilla JS + HTML + CSS in the browser, Node.js +
Express API, PostgreSQL, deployable to AWS. 163 tests in ~1.5 seconds.

**A course.** [`docs/`](docs/00-start-here.md) contains twelve chapters
explaining why every part of it is shaped the way it is — architecture, OOP,
domain modelling, SQL, REST, HTML/CSS, the DOM, testing, maintenance and
deployment. The comments in the source are part of the teaching material.

→ **[Start the course](docs/00-start-here.md)**

## No build step

The browser loads ES modules natively. The server serves them as files. There is
no bundler, no transpiler, and no config file doing something invisible on your
behalf — **the code that runs is the code you wrote.**

Five dependencies in total: `express`, `cors`, `pg`, `dotenv`, and `jsdom` for
tests.

## Run it

```bash
npm install
npm run db:up      # PostgreSQL in Docker
npm run db:migrate # create the tables
npm run db:seed    # 124 cited events across 12 countries
npm run dev        # http://localhost:4000
```

One server, one URL — it serves both the API and the web page, so there is no
proxy and no CORS to configure while you are learning.

No Docker? Install PostgreSQL locally:

```bash
createdb history_timeline
echo 'DATABASE_URL=postgres://localhost:5432/history_timeline' > .env
npm run db:migrate && npm run db:seed && npm run dev
```

```bash
npm test                # 163 tests
npm run test:coverage   # where the blind spots are
```

## The MVP features

1. **A timeline** of years, 1900 → present, with bar height showing how much is
   recorded in each year
2. **Select a year** — by click, or arrow / Page Up / Page Down / Home / End
3. **List the countries** that have events in that year, grouped by region
4. **Read the events** for that country in that year, each with its source

The domain already supports **100 AD onwards**; the MVP only *renders* from 1900.
Widening it is a change to a default, not to the architecture.

## Architecture

```
views → App (state) → TimelineApi → HTTP
                                     ↓
      controllers → TimelineService → Repository (abstract)
                                        ├── PostgresCountryRepository
                                        └── InMemoryCountryRepository (tests)
                                             ↓
                                        domain: Year, Country, HistoricalEvent
```

Dependencies point inward. The domain imports nothing. The service knows neither
Express nor `pg` — which is why the server suite runs in under a second against
in-memory repositories, and why the same code runs on PostgreSQL in production.

In the browser: **views never call `fetch`**, exactly as services never write SQL.

## Layout

```
docs/          the course, chapters 00–11
server/
  src/domain/       rules: Year, Country, HistoricalEvent
  src/repositories/ data access — abstract + two implementations
  src/services/     use cases
  src/api/          HTTP: controllers, error mapping
  src/db/           connection, migrations, seed data
  src/main.js       composition root
client/          served as static files; no build
  index.html        the page structure, written by hand
  css/styles.css    one stylesheet, custom properties, dark mode
  js/dom.js         30-line element builder (and the XSS defence)
  js/Store.js       40-line observable state container
  js/api/           ApiClient base class, TimelineApi
  js/components/    TimelineView, CountryListView, EventListView
  js/App.js         owns state, loads data
infra/         AWS deployment
```

## API

| Method & path | Returns |
|---|---|
| `GET /health` · `GET /ready` | liveness · readiness (checks the database) |
| `GET /api/timeline?startYear=&endYear=` | `{ startYear, endYear, ticks[] }` |
| `GET /api/countries` | `{ countries[] }` |
| `GET /api/years/:year/countries` | `{ year, countries[] }` |
| `GET /api/years/:year/events` | `{ year, events[] }` |
| `GET /api/years/:year/countries/:code/events` | `{ year, country, events[] }` |

```bash
curl localhost:4000/api/years/1994/countries/ZAF/events
```

## The editorial rules

This is a factual history project, and the rules are enforced in code:

1. **Factual entries only.**
2. **Every event cites a source** — `HistoricalEvent` refuses to be constructed
   without a valid URL, so an uncited event cannot exist anywhere in the system,
   not even in a test.
3. **Neutral summaries** — describe what happened; leave judgement to the reader.

Rule 2 is the best illustration in the project of what a domain model is *for*.

## Contributing an event

Add it to `server/src/db/seedData.js`, then `npm run db:seed` (idempotent — safe
to re-run). `npm test` will tell you if it breaks an editorial rule: a missing
source, a duplicate, an unknown country code, a year outside the window.

## Licence

The code is yours to use. Event summaries are original prose; each links to the
source that supports it.
