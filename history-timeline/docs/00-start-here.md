# Start here

This repository is a **finished, working project that is also a course**. The
code in `server/` and `client/` runs; the chapters in `docs/` explain why every
part of it is shaped the way it is.

## What you will build

A website with a timeline. You pick a year, you pick a country, and you read
what is documented to have happened there — with a link to the source for every
claim.

The four MVP features:

| # | Feature | Where it lives |
|---|---------|----------------|
| 1 | A timeline showing years (1900 → present) | `Timeline.tsx` → `GET /api/timeline` → `TimelineService.getTimelineScale` |
| 2 | Select a specific year | `App.tsx` state → every panel derives from it |
| 3 | List the countries | `CountryPicker.tsx` → `GET /api/years/:year/countries` |
| 4 | Read that country's events for that year | `EventList.tsx` → `GET /api/years/:year/countries/:code/events` |

The domain already supports **100 AD onwards** (`Year.EARLIEST = 100`). The MVP
only *renders* from 1900 (`Year.MVP_START`). Widening the window later is a
change to a default, not to the architecture — that is deliberate, and chapter 3
explains why designing for it now cost nothing.

## How to use this course

Each chapter follows the same rhythm:

> **Concept** → **The code in this project** → **Do it yourself** → **Check yourself**

Read the chapter, then open the file it names and read the real code with its
comments. The comments in the source are part of the teaching material, not
decoration. Then do the exercise at the end — that is where the learning
actually happens.

### Two ways through

- **Follow along** — read the chapters in order and study the working code.
  Fastest route to understanding the whole stack.
- **Build it yourself** — delete `server/src` and `client/src`, keep `docs/` and
  the tests, and make the tests pass one chapter at a time. Much slower, much
  more durable. The tests are written to be a specification.

## The chapters

| Chapter | Topic | You will learn |
|---------|-------|----------------|
| [01](01-architecture.md) | Architecture | Layers, dependency direction, why this shape |
| [02](02-typescript-and-oop.md) | TypeScript & OOP | Classes, interfaces, SOLID as used here |
| [03](03-domain-model.md) | The domain model | Entities vs value objects, making bad states impossible |
| [04](04-postgres.md) | PostgreSQL | Schema design, indexes, migrations, SQL injection |
| [05](05-repositories.md) | The repository pattern | Dependency inversion, contract tests |
| [06](06-http-api.md) | The HTTP API | REST design, status codes, error handling |
| [07](07-react-frontend.md) | React | Components, hooks, state design, race conditions |
| [08](08-testing.md) | Testing | The pyramid, what to test, what never to test |
| [09](09-maintenance.md) | Maintaining it | CI, logging, dependencies, code review |
| [10](10-aws-deployment.md) | AWS | Deploying the whole thing, and what it costs |

## Prerequisites

- **Node.js 20+** and npm — `node -v`
- **Docker** (for PostgreSQL) or a local PostgreSQL 14+ install
- Comfort with JavaScript basics: variables, functions, `async`/`await`, arrays.
  You do **not** need prior TypeScript, React, SQL or AWS experience.

## Run it now

Before reading anything, get it working. Five commands:

```bash
cd history-timeline
npm install                 # installs both workspaces
npm run db:up               # starts PostgreSQL in Docker
npm run db:migrate          # creates the tables
npm run db:seed             # loads 124 cited events across 12 countries
npm run dev                 # API on :4000, site on http://localhost:5173
```

Open http://localhost:5173, click a tall bar on the timeline, pick a country.

**No Docker?** Install PostgreSQL locally, create the database, and point
`server/.env` at it:

```bash
createdb history_timeline
echo 'DATABASE_URL=postgres://localhost:5432/history_timeline' > server/.env
npm run db:migrate && npm run db:seed && npm run dev
```

Run the tests too — they should be green before you change anything:

```bash
npm test        # 132 tests, about two seconds
```

## The shape of the repository

```
history-timeline/
├── docs/               ← the course you are reading
├── server/             ← Node.js + Express + PostgreSQL API
│   ├── src/
│   │   ├── domain/         the rules: Year, Country, HistoricalEvent
│   │   ├── repositories/   data access, one abstract class per aggregate
│   │   ├── services/       use cases (TimelineService)
│   │   ├── api/            HTTP: controllers, error handling
│   │   ├── db/             connection, migrations, seed data
│   │   ├── config/         validated environment
│   │   └── main.ts         the composition root
│   └── tests/          unit + integration
├── client/             ← React + TypeScript + Vite
│   ├── src/
│   │   ├── api/            ApiClient base class, TimelineApi
│   │   ├── hooks/          useAsyncResource and friends
│   │   ├── components/     Timeline, CountryPicker, EventList
│   │   └── App.tsx         state owner
│   └── tests/
├── infra/              ← AWS deployment
└── docker-compose.yml  ← local PostgreSQL
```

## A word about the subject matter

This is a **factual history** project, and that constrains the engineering. The
editorial rules live in `server/src/db/seedData.ts` and are enforced in code:

1. Factual entries only.
2. **Every event cites a source.** `HistoricalEvent` refuses to be constructed
   without a valid URL — an uncited event cannot exist in this system, not even
   in a test.
3. Neutral summaries: describe what happened, leave judgement to the reader.

Rule 2 is the best illustration in the project of what a domain model is *for*.
It is not a data-shaping convenience. It is where the thing your product
promises gets turned into something the compiler and the test suite enforce.

→ Next: [01 — Architecture](01-architecture.md)
