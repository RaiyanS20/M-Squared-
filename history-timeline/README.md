# History Timeline

A factual, educational history site — **and a project-based course in building
one**.

Pick a year on the timeline, choose a country, and read what is documented to
have happened there. Every event links to its source.

![Stack](https://img.shields.io/badge/React-19-61dafb) ![Node](https://img.shields.io/badge/Node-22-339933) ![TypeScript](https://img.shields.io/badge/TypeScript-5.9-3178c6) ![Postgres](https://img.shields.io/badge/PostgreSQL-16-336791)

---

## Two things at once

**A working application.** React + TypeScript frontend, Node.js + Express API,
PostgreSQL, deployable to AWS. 132 tests, ~3 seconds.

**A course.** [`docs/`](docs/00-start-here.md) contains eleven chapters that
explain why every part of it is shaped the way it is — architecture, OOP, domain
modelling, SQL, REST, React, testing, maintenance and deployment. The comments in
the source code are part of the teaching material.

→ **[Start the course](docs/00-start-here.md)**

## Run it

```bash
npm install        # both workspaces
npm run db:up      # PostgreSQL in Docker
npm run db:migrate # create the tables
npm run db:seed    # 124 cited events across 12 countries
npm run dev        # API :4000, site http://localhost:5173
```

No Docker? Install PostgreSQL locally, then:

```bash
createdb history_timeline
echo 'DATABASE_URL=postgres://localhost:5432/history_timeline' > server/.env
npm run db:migrate && npm run db:seed && npm run dev
```

```bash
npm test           # 132 tests
npm run build      # compile both workspaces
```

## The MVP features

1. **A timeline** of years, 1900 → present, with bar height showing how much is
   recorded in each year
2. **Select a year** — by click or by keyboard
3. **List the countries** that have events in that year, grouped by region
4. **Read the events** for that country in that year, each with its source

The domain already supports **100 AD onwards**; the MVP only *renders* from 1900.
Widening it is a change to a default, not to the architecture.

## Architecture

```
React components → hooks → TimelineApi → HTTP
                                          ↓
  controllers → TimelineService → Repository (abstract)
                                     ├── PostgresCountryRepository
                                     └── InMemoryCountryRepository (tests)
                                          ↓
                                     domain: Year, Country, HistoricalEvent
```

Dependencies point inward. The domain imports nothing. The service knows neither
Express nor `pg` — which is why the whole server suite runs in about a second
against in-memory repositories, and why the same code runs on PostgreSQL in
production.

## Layout

```
docs/          the course, chapters 00–10
server/        Node.js + Express + PostgreSQL
  src/domain/       rules: Year, Country, HistoricalEvent
  src/repositories/ data access — abstract + two implementations
  src/services/     use cases
  src/api/          HTTP: controllers, error mapping
  src/db/           connection, migrations, seed data
  src/main.ts       composition root
client/        React + TypeScript + Vite
  src/api/          ApiClient base class, TimelineApi
  src/hooks/        useAsyncResource and friends
  src/components/   Timeline, CountryPicker, EventList
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

Add it to `server/src/db/seedData.ts`, then `npm run db:seed` (idempotent — safe
to re-run). `npm test` will tell you if it breaks an editorial rule: a missing
source, a duplicate, an unknown country code, a year outside the window.

## Licence

The code is yours to use. Event summaries are original prose; each links to the
source that supports it.
