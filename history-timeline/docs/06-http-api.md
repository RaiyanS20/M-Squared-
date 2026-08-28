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

Plus everything else: the static browser app.

## REST design decisions

**Paths name resources; they nest to show containment.**
`/api/years/1994/countries/ZAF/events` reads as a sentence: the events, of South
Africa, in 1994. Compare `/api/getEvents?y=1994&c=ZAF` — a function call wearing
a URL costume.

**Codes in the URL, not ids.** `/1969/FRA` beats `/1969/74`: readable, stable,
shareable. Ids are an implementation detail of our database; ISO 3166 codes are a
fact about the world.

**Objects at the top level, not bare arrays.** `{ countries: [...] }`, never
`[...]`. Adding `{ countries, totalCount }` later is then a *non-breaking*
change. Returning a bare array paints you into a corner on day one.

**Nothing found is `200 []`, not `404`.** Asking for France in 1955 is a valid
question with the answer "nothing recorded". `404` means *the resource does not
exist* — which is what `/countries/ZZZ` gets, because there is no such country.
This distinction confuses people constantly; the test suite pins both cases.

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

```js
#listEventsForYearAndCountry = async (req, res) => {
  const year = this.#requiredYear(req.params.year);          // 1. read + validate
  const code = String(req.params.code ?? '');
  if (!/^[A-Za-z]{3}$/.test(code)) throw new ValidationError(...);
  res.json(await this.#service.getEventsForYearAndCountry(year, code));  // 2 + 3
};
```

**The moment a controller contains an `if` about history, that rule is in the
wrong file.** Validating *shape* ("is this three letters?") belongs here;
validating *meaning* ("does this country exist?") belongs in the service.

### Why arrow-function fields

```js
#getTimeline = async (req, res) => { ... };   // arrow field, not a method
```

A plain method loses its `this` when detached:

```js
router.get('/timeline', this.getTimeline);
// as a method: `this` is undefined inside when Express calls it → TypeError
```

Arrow fields capture `this` at construction, so the reference is safe to pass
around. Chapter 2 covers this in full; it is the most common `this` bug in
JavaScript.

### A small detail worth copying

```js
#optionalYear(raw, field) {
  if (raw === undefined || raw === '') return undefined;
  if (typeof raw !== 'string') throw new ValidationError(`${field} must be a single value.`, field);
  return Year.fromString(raw).value;
}
```

`?startYear=1&startYear=2` makes Express hand you an **array**. Without that
check it becomes a confusing failure deeper in the stack. There is a test for it.

## Error handling in one place

```js
if (error instanceof DomainError) {
  res.status(error.httpStatus).json({ error: { code: error.code, message: error.message } });
  return;
}
console.error('[api] unhandled error', error);
res.status(500).json({ error: { code: 'INTERNAL_ERROR',
  message: isProduction ? 'Something went wrong.' : ... } });
```

Three things worth copying:

1. **Polymorphic mapping.** The handler never names `NotFoundError`. Add a new
   `DomainError` subclass and it is handled — open/closed in production code.
2. **Anything not a `DomainError` is a bug.** Log it fully, tell the client
   nothing. A stack trace tells an attacker your file layout, library versions
   and often your queries.
3. **A machine-readable `code` beside the human `message`.** The client switches
   on `NOT_FOUND`, never on English text — so messages can be reworded or
   translated freely.

**Express identifies error middleware by arity: a function of four arguments.**
Write three and it is registered as ordinary middleware that never sees an error.
That is a genuinely baffling bug the first time you hit it.

## `createApp` vs `main.js`

```js
export function createApp(deps) { ... }   // builds
app.listen(config.port, ...)               // runs — main.js only
```

**Separating "build the app" from "run the app" is the single change that makes
an HTTP layer testable.** The tests call `createApp` with in-memory
repositories:

```js
const app = createApp({ service, db: null, isProduction: false, serveClient: false });
server = app.listen(0, resolve);   // port 0 = "any free port"
```

Note `serveClient: false` in tests: with the static handler on, an unknown path
would fall through to it instead of returning the JSON 404 the test is checking.

## Middleware order

Order is behaviour, not style:

```js
app.use(cors(...));               // 1. before routes
app.use(express.json(...));       // 2. parse bodies before handlers read them
app.use('/', health.routes());    // 3. routes
app.use('/api', timeline.routes());
app.use(express.static(CLIENT_DIR));   // 4. the browser app
app.use(notFoundHandler);         // 5. nothing matched
app.use(errorHandler(...));       // 6. LAST
```

Put the error handler before the routes and it never fires. Put
`notFoundHandler` before the routes and *every* request is a 404. Put
`express.static` before `/api` and a file named `api` would shadow your API.

## CORS, honestly

The browser refuses cross-origin requests unless the server opts in.

**In this project, development has no CORS at all** — the API serves the browser
app, so everything is one origin. The configuration exists for production, where
the site is on CloudFront and the API is behind a load balancer:

```js
if (corsOrigins.length > 0) {
  app.use(cors({ origin: corsOrigins, methods: ['GET'] }));
}
```

**Never `origin: '*'`.** That tells every website on the internet that their
JavaScript may read your responses. An empty list denies all, which is the right
default.

## Health vs readiness

```js
GET /health  → { status: 'ok' }        // is the process alive?
GET /ready   → checks the database, 503 if down
```

`/health` answers "should this be restarted?". `/ready` answers "should this
receive traffic?" — a container with a dead database connection is pulled from
the load balancer instead of serving 500s. Point the AWS health check at
`/ready`.

`HealthController` takes anything with an `isHealthy()` method, not a
`PostgresDatabase`. That is interface segregation, and it keeps the API layer
free of any storage type — the tests pass `null`.

## Small security details

```js
app.set('trust proxy', true);      // read X-Forwarded-For from the load balancer
app.disable('x-powered-by');       // stop advertising "Express" to scanners
express.json({ limit: '100kb' })   // a 2GB body should be rejected, not buffered
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
  exist. Clients handle the two differently — one shows *"nothing yet"*, the
  other *"that country isn't in our data"*.
- **`code` plus `message`**: `code` is a stable contract for the client to branch
  on; `message` is for humans and free to change. Branching on message text
  breaks the moment someone fixes a typo.
- **Four arguments**: Express identifies error middleware by `fn.length === 4`
  (`err, req, res, next`). Write three and it never sees an error.

</details>

→ Next: [07 — HTML and CSS](07-html-and-css.md)
