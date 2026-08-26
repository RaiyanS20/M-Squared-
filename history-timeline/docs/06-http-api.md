# 06 — The HTTP API

> **Concept** → The code → Do it yourself → Check yourself

## The routes

| Method & path | Feature | Returns |
|---|---|---|
| `GET /health` | liveness | `{ status, uptimeSeconds }` |
| `GET /ready` | readiness | `{ status, database }` — 503 if the DB is down |
| `GET /api/timeline?startYear=&endYear=` | 1 | `{ startYear, endYear, ticks[] }` |
| `GET /api/countries` | 3 | `{ countries[] }` |
| `GET /api/years/:year/countries` | 2+3 | `{ year, countries[] }` |
| `GET /api/years/:year/events` | — | `{ year, events[] }` |
| `GET /api/years/:year/countries/:code/events` | 4 | `{ year, country, events[] }` |

## REST design decisions

**Paths name resources; they nest to show containment.**
`/api/years/1994/countries/ZAF/events` reads as a sentence: the events, of South
Africa, in 1994. Compare `/api/getEvents?y=1994&c=ZAF` — that is a function call
wearing a URL costume.

**Codes in the URL, not ids.** `/1969/FRA` beats `/1969/74`: readable, stable,
and shareable. Ids are an implementation detail of our database; ISO 3166 codes
are a fact about the world.

**Objects at the top level, not bare arrays.** `{ countries: [...] }`, never
`[...]`. Adding `{ countries, totalCount }` later is then a *non-breaking*
change. Returning a bare array paints you into a corner on day one.

**Nothing found is `200 []`, not `404`.** Asking for France in 1955 is a
perfectly valid question with the answer "nothing recorded". `404` means *the
resource does not exist* — which is what `/countries/ZZZ` gets, because there is
no such country. This distinction confuses people constantly; the test suite
pins both cases.

## Status codes that matter here

| Code | When | Example |
|---|---|---|
| `200` | Fine, including empty results | France in 1955 |
| `400` | The request itself is malformed | `?startYear=abc`, year 3000, code `FRANCE` |
| `404` | Resource does not exist | country `ZZZ`, unknown route |
| `500` | We have a bug | anything unanticipated |
| `503` | Alive but not serving | `/ready` with the database down |

The rule of thumb: **4xx means the caller must change something; 5xx means we
must.** If a client can retry unchanged and succeed, it was not a 4xx.

## Thin controllers

Every controller method does exactly three things: read input, call one service
method, shape the response.

```ts
private listEventsForYearAndCountry = async (req, res) => {
  const year = this.requiredYear(req.params.year);   // 1. read + validate
  const code = String(req.params.code ?? '');
  if (!/^[A-Za-z]{3}$/.test(code)) throw new ValidationError(...);
  const result = await this.service.getEventsForYearAndCountry(year, code);  // 2.
  res.json(result);                                                          // 3.
};
```

**The moment a controller contains an `if` about history, that rule is in the
wrong file.** Validation of *shape* ("is this three letters?") belongs here;
validation of *meaning* ("does this country exist?") belongs in the service.

### Why arrow-function properties

```ts
private getTimeline = async (req, res) => { ... };   // arrow property
```

Not a method. A plain method loses its `this` when detached:

```ts
router.get('/timeline', this.getTimeline);
// as a method: `this` inside is undefined when Express calls it → TypeError
```

Arrow properties capture `this` at construction, so the reference is safe to pass
around. (`.bind(this)` in the constructor is the equivalent older fix.) This is
one of the most common "why is `this` undefined" bugs in TypeScript, and it is
worth understanding rather than pattern-matching.

## Error handling in one place

```ts
if (error instanceof DomainError) {
  res.status(error.httpStatus).json({ error: { code: error.code, message: error.message } });
  return;
}
console.error('[api] unhandled error', error);
res.status(500).json({ error: { code: 'INTERNAL_ERROR', message: isProduction ? 'Something went wrong.' : ... } });
```

Three things worth copying:

1. **Polymorphic mapping.** The handler never names `NotFoundError`. Add a new
   `DomainError` subclass and it is handled — that is open/closed in production
   code rather than in a textbook.
2. **Anything not a `DomainError` is a bug.** Log it fully, tell the client
   nothing. In production a stack trace tells an attacker your file layout,
   library versions and often your queries.
3. **A machine-readable `code` beside the human `message`.** The client switches
   on `NOT_FOUND`, never on English text. This is why `ApiError` on the frontend
   carries `code` — and why messages can be reworded or translated freely.

The shape is consistent for *every* error, including unknown routes — so the
React app never has to handle an HTML error page.

## `createApp` vs `main.ts`

```ts
export function createApp(deps: AppDependencies): Express   // builds
app.listen(config.port, ...)                                 // runs — main.ts only
```

**Separating "build the app" from "run the app" is the single change that makes
an HTTP layer testable.** Tests call `createApp` with in-memory repositories and
drive it with supertest — no port, no database, no cleanup:

```ts
const app = createApp({ service, db: null, corsOrigins: [...], isProduction: false });
await request(app).get('/api/years/1969/countries/FRA/events').expect(200);
```

22 integration tests, no infrastructure.

## Middleware order

Order is not stylistic; it is behaviour:

```ts
app.use(cors(...));            // 1. before routes — must run for every request
app.use(express.json(...));    // 2. parse bodies before handlers read them
app.use('/', health.routes()); // 3. routes
app.use('/api', timeline.routes());
app.use(notFoundHandler);      // 4. nothing matched
app.use(errorHandler(...));    // 5. LAST — four args marks it as the error handler
```

Put the error handler before the routes and it never fires. Put `notFoundHandler`
before the routes and *every* request is a 404.

## CORS, honestly

The browser refuses cross-origin requests unless the server opts in. In
development the site is `localhost:5173` and the API is `localhost:4000` —
different origins.

```ts
cors({ origin: deps.corsOrigins.length > 0 ? deps.corsOrigins : false, methods: ['GET'] })
```

**Never `origin: '*'` on an API with credentials.** It tells every website on the
internet that their JavaScript may read your responses. The allow-list comes from
`CORS_ORIGINS`, and there is a test asserting an unknown origin is refused.

Note also that the Vite dev server proxies `/api` to `:4000` (see
`vite.config.ts`), so in day-to-day development the browser sees one origin and
CORS never triggers. The configuration matters in production, where the site and
API really are on different hosts.

## Health vs readiness

```ts
GET /health  → { status: 'ok' }              // is the process alive?
GET /ready   → checks the database, 503 if down
```

The distinction matters on AWS. `/health` is cheap and answers "should this task
be restarted?". `/ready` answers "should this task receive traffic?" — a task
with a dead database connection is pulled from the load balancer instead of
serving 500s. Point the ALB health check at `/ready`.

## Small security details in `createApp`

```ts
app.set('trust proxy', true);   // read X-Forwarded-For from the ALB, so req.ip is the real client
app.disable('x-powered-by');    // stop advertising "Express" to scanners
express.json({ limit: '100kb' }) // a 2GB body should be rejected, not buffered
```

## Do it yourself

1. **Meet the errors.** With the server running:
   ```bash
   curl -i localhost:4000/api/years/1969/countries/ZZZ/events   # 404 NOT_FOUND
   curl -i localhost:4000/api/years/3000/countries              # 400 VALIDATION_ERROR
   curl -i localhost:4000/api/years/1955/countries/FRA/events   # 200, events: []
   curl -i localhost:4000/api/nope                              # 404 ROUTE_NOT_FOUND
   ```
   Explain why the third is not a 404.

2. **Add an endpoint.** `GET /api/countries/:code` returning one country, 404 if
   unknown. Service method, controller route, integration test. Notice you do not
   have to touch the error handler.

3. **Break the middleware order.** Move `app.use(errorHandler(...))` above the
   routes. Watch `curl localhost:4000/api/years/3000/countries` return an HTML
   stack trace instead of JSON. Put it back.

4. **Add pagination.** `GET /api/years/:year/events?limit=&offset=`. Where does
   `limit` get validated — controller or service? (Both, arguably: shape in the
   controller, a maximum in the service. Decide and defend it.)

## Check yourself

- Why is `/api/years/1955/countries/FRA/events` a 200 and not a 404?
- Why does the error response carry `code` as well as `message`?
- Why must `errorHandler` take four arguments?

<details>
<summary>Answers</summary>

- **200 with an empty list**: the resource (France, 1955) exists and the question
  is valid; the answer is "nothing recorded". 404 would mean France does not
  exist. Clients handle "empty" and "missing" differently — one shows *"nothing
  yet"*, the other shows *"that country isn't in our data"*.
- **`code` plus `message`**: `code` is a stable contract for the client to branch
  on; `message` is for humans and free to change. Branching on message text
  breaks the moment someone fixes a typo.
- **Four arguments**: Express identifies error-handling middleware by
  `fn.length === 4` (`err, req, res, next`). Write three and it is registered as
  ordinary middleware that never sees an error.

</details>

→ Next: [07 — The React frontend](07-react-frontend.md)
