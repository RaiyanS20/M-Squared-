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
| 1 | A timeline showing years (1900 → present) | `TimelineView.js` → `GET /api/timeline` → `TimelineService.getTimelineScale` |
| 2 | Select a specific year | `App.js` state → every panel derives from it |
| 3 | List the countries | `CountryListView.js` → `GET /api/years/:year/countries` |
| 4 | Read that country's events for that year | `EventListView.js` → `GET /api/years/:year/countries/:code/events` |

The domain already supports **100 AD onwards** (`Year.EARLIEST = 100`). The MVP
only *renders* from 1900 (`Year.MVP_START`). Widening the window later is a
change to a default, not to the architecture.

## The stack, and why

| Layer | Choice | Why |
|---|---|---|
| Browser | **Vanilla JavaScript, HTML, CSS** | Learn the platform before a framework. Everything React does for you, you will do once by hand — and then you will know what you are buying. |
| Server | **Node.js + Express** | Same language on both sides, so there is one language to learn, not two. Express is the Node standard. |
| Database | **PostgreSQL** | The data is relational, and Postgres enforces rules the application cannot bypass. |
| Cloud | **AWS** | Covered in chapter 11, with honest costs. |
| Style | **OOP with ES2022 classes** | Real `#private` fields, `extends`, static factories — as used in industry. |

**No build step anywhere.** The browser loads ES modules natively; the server
serves them as files. No bundler, no transpiler, no config file doing something
invisible on your behalf. When something breaks, the thing that broke is code
you wrote.

## How to use this course

Each chapter follows the same rhythm:

> **Concept** → **The code in this project** → **Do it yourself** → **Check yourself**

Read the chapter, then open the file it names and read the real code with its
comments. The comments in the source are part of the teaching material. Then do
the exercise — that is where the learning actually happens.

### Two ways through

- **Follow along** — read in order and study the working code. Fastest route to
  understanding the whole stack.
- **Build it yourself** — delete `server/src` and `client/js`, keep `docs/` and
  the tests, and make the tests pass one chapter at a time. Much slower, much
  more durable. The tests are written to be a specification.

## The chapters

| Chapter | Topic | You will learn |
|---------|-------|----------------|
| [01](01-architecture.md) | Architecture | Layers, dependency direction, why this shape |
| [02](02-javascript-and-oop.md) | JavaScript & OOP | Modules, classes, `#private`, `this`, SOLID |
| [03](03-domain-model.md) | The domain model | Entities vs value objects, validation as a type system |
| [04](04-postgres.md) | PostgreSQL | Schema design, indexes, migrations, SQL injection |
| [05](05-repositories.md) | The repository pattern | Dependency inversion, contract tests |
| [06](06-http-api.md) | The HTTP API | REST design, status codes, error handling |
| [07](07-html-and-css.md) | HTML & CSS | Semantic markup, accessibility, the cascade, layout |
| [08](08-javascript-in-the-browser.md) | JS in the browser | The DOM, XSS, state, events, race conditions |
| [09](09-testing.md) | Testing | The pyramid, what to test, what never to test |
| [10](10-maintenance.md) | Maintaining it | CI, logging, dependencies, code review |
| [11](11-aws-deployment.md) | AWS | Deploying the whole thing, and what it costs |

## Prerequisites

- **Node.js 20+** — check with `node -v`
- **Docker** (for PostgreSQL) or a local PostgreSQL 14+ install
- Basic programming familiarity: variables, functions, loops, arrays. You do
  **not** need prior JavaScript, SQL, or AWS experience.

## Run it now

Before reading anything, get it working:

```bash
cd history-timeline
npm install                 # four runtime dependencies
npm run db:up               # starts PostgreSQL in Docker
npm run db:migrate          # creates the tables
npm run db:seed             # loads 124 cited events across 12 countries
npm run dev                 # http://localhost:4000
```

Open <http://localhost:4000>, click a tall bar on the timeline, pick a country.

There is **one** server and **one** URL — it serves both the API and the web
page. That is a deliberate simplification: no separate dev server, no proxy, no
CORS to configure while you are learning everything else.

**No Docker?** Install PostgreSQL locally:

```bash
createdb history_timeline
echo 'DATABASE_URL=postgres://localhost:5432/history_timeline' > .env
npm run db:migrate && npm run db:seed && npm run dev
```

Run the tests too — they should be green before you change anything:

```bash
npm test        # 163 tests, about 1.5 seconds
```

## The shape of the repository

```
history-timeline/
├── docs/                ← the course you are reading
├── server/
│   ├── src/
│   │   ├── domain/          the rules: Year, Country, HistoricalEvent
│   │   ├── repositories/    data access, one abstract base per aggregate
│   │   ├── services/        use cases (TimelineService)
│   │   ├── api/             HTTP: controllers, error handling
│   │   ├── db/              connection, migrations, seed data
│   │   ├── config/          validated environment
│   │   └── main.js          the composition root
│   └── tests/
├── client/              ← no build step; served as static files
│   ├── index.html           the page structure, written by hand
│   ├── css/styles.css       one stylesheet, custom properties
│   ├── js/
│   │   ├── dom.js             a 30-line element builder
│   │   ├── Store.js           a 40-line observable state container
│   │   ├── api/               ApiClient base class, TimelineApi
│   │   ├── components/        TimelineView, CountryListView, EventListView
│   │   ├── App.js             owns state, loads data
│   │   └── app.js             entry point
│   └── tests/
├── infra/               ← AWS deployment
└── docker-compose.yml   ← local PostgreSQL
```

## A word about the subject matter

This is a **factual history** project, and that constrains the engineering. The
editorial rules live in `server/src/db/seedData.js` and are enforced in code:

1. Factual entries only.
2. **Every event cites a source.** `HistoricalEvent` refuses to be constructed
   without a valid URL — an uncited event cannot exist in this system, not even
   in a test.
3. Neutral summaries: describe what happened, leave judgement to the reader.

Rule 2 is the best illustration in the project of what a domain model is *for*.
It is not a data-shaping convenience. It is where the thing your product
promises becomes something the test suite enforces.

→ Next: [01 — Architecture](01-architecture.md)
